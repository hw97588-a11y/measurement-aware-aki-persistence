// Fast, fixed-reference controlled thinning simulation. See the paired
// Python cache/aggregation scripts for definitions and non-causal constraints.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

constexpr double M48 = 48.0 * 60.0;
constexpr double M7D = 7.0 * 24.0 * 60.0;
struct Point { double t; double c; };
struct Episode { int hospital; int category; double end; std::vector<Point> values; };

std::vector<std::string> split(const std::string& text, char sep) {
    std::vector<std::string> result; std::string part; std::stringstream stream(text);
    while (std::getline(stream, part, sep)) result.push_back(part);
    return result;
}

std::vector<Episode> read_cache(const std::string& file, int& hospitals) {
    std::ifstream input(file); if (!input) throw std::runtime_error("cannot open cache");
    std::vector<Episode> episodes; std::string line; hospitals = 0;
    while (std::getline(input, line)) {
        auto fields = split(line, '\t'); if (fields.size() != 4) continue;
        Episode e; e.hospital = std::stoi(fields[0]); e.category = std::stoi(fields[1]); e.end = std::stod(fields[2]); hospitals = std::max(hospitals, e.hospital + 1);
        for (const auto& token : split(fields[3], ';')) { auto pair = split(token, ','); if (pair.size() == 2) e.values.push_back({std::stod(pair[0]), std::stod(pair[1])}); }
        episodes.push_back(std::move(e));
    }
    return episodes;
}

std::vector<Point> thin(const std::vector<Point>& points, double phase, double interval, int rule) {
    std::vector<Point> output; long current = 0; bool have = false; double best = 0.0;
    for (const auto& point : points) {
        long bin = static_cast<long>(std::floor((point.t - phase) / interval));
        double centre = phase + (static_cast<double>(bin) + .5) * interval;
        double distance = std::abs(point.t - centre);
        if (!have || bin != current) { output.push_back(point); current = bin; best = distance; have = true; }
        else if (rule == 2 || (rule == 0 && (distance < best || (distance == best && point.t < output.back().t)))) { output.back() = point; best = distance; }
    }
    return output;
}

// -1 not detected; 0 transient; 1 persistent; 2 interval; 3 right-censored;
// 4 detected but no longer has 48-h potential coverage.
int classify(const std::vector<Point>& s, double end) {
    bool previous_positive = false; bool have_last_non = false; double last_non = 0.0;
    for (size_t i = 0; i < s.size(); ++i) {
        double time = s[i].t, value = s[i].c; bool have48 = false, have7 = false; double base48 = 1e100, base7 = 1e100;
        for (size_t k = 0; k < i; ++k) { double gap = time - s[k].t; if (gap > 0 && gap <= M7D) { have7 = true; base7 = std::min(base7, s[k].c); if (gap <= M48) { have48 = true; base48 = std::min(base48, s[k].c); } } }
        bool positive = (have48 && value - base48 >= .3) || (have7 && value / base7 >= 1.5);
        if (positive && !previous_positive && time >= 0.0 && time <= std::min(M7D, end) && have_last_non && have7) {
            double recovery_limit = std::min(base7 + .3, 1.5 * base7); int recovery = -1;
            for (size_t j = i; j < s.size(); ++j) if (s[j].c < recovery_limit) { recovery = static_cast<int>(j); break; }
            int category;
            if (recovery >= 0) {
                double lower = std::max(0.0, s[recovery - 1 < static_cast<int>(i) ? i : recovery - 1].t - time);
                double upper = s[recovery].t - last_non;
                category = upper <= M48 ? 0 : (lower > M48 ? 1 : 2);
            } else {
                double last_positive = time; for (size_t j = i; j < s.size(); ++j) if (s[j].c >= recovery_limit) last_positive = s[j].t;
                category = (last_positive - time > M48) ? 1 : 3;
            }
            return end >= time + M48 ? category : 4;
        }
        if (!positive) { have_last_non = true; last_non = time; }
        previous_positive = positive;
    }
    return -1;
}

std::vector<double> ranks(const std::vector<double>& values) {
    std::vector<int> order(values.size()); std::iota(order.begin(), order.end(), 0); std::sort(order.begin(), order.end(), [&](int a, int b){ return values[a] < values[b]; });
    std::vector<double> output(values.size()); size_t start = 0;
    while (start < order.size()) { size_t stop = start + 1; while (stop < order.size() && values[order[stop]] == values[order[start]]) ++stop; double rank = (static_cast<double>(start + 1) + static_cast<double>(stop)) / 2.0; for (size_t j = start; j < stop; ++j) output[order[j]] = rank; start = stop; }
    return output;
}
double correlation(const std::vector<double>& a, const std::vector<double>& b) {
    double ma = std::accumulate(a.begin(), a.end(), 0.0) / a.size(), mb = std::accumulate(b.begin(), b.end(), 0.0) / b.size(), xy=0, xa=0, xb=0;
    for (size_t i=0;i<a.size();++i) { double da=a[i]-ma, db=b[i]-mb; xy+=da*db; xa+=da*da; xb+=db*db; }
    return xy / std::sqrt(xa*xb);
}

int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "usage: sim cache output_dir replicates seed\n"; return 2; }
    int hospitals=0, replicates=std::stoi(argv[3]); auto episodes = read_cache(argv[1], hospitals); std::string outdir=argv[2]; std::mt19937_64 rng(std::stoull(argv[4]));
    std::vector<int> hcount(hospitals,0), hraw(hospitals,0); int refcat[4]{};
    for (const auto& e: episodes) { ++hcount[e.hospital]; hraw[e.hospital] += e.category >= 2; ++refcat[e.category]; }
    std::vector<int> selected; std::vector<double> rawrate; for (int h=0;h<hospitals;++h) if (hcount[h]>=20) { selected.push_back(h); rawrate.push_back(static_cast<double>(hraw[h])/hcount[h]); }
    auto rawrank=ranks(rawrate);
    for (int rule=0; rule<3; ++rule) for (int hours: {12,24,36,48}) {
        std::string name = rule==0 ? "nearest" : (rule==1 ? "first" : "last");
        std::ofstream metrics(outdir + "/metrics_" + name + "_" + std::to_string(hours) + ".tsv"); std::ofstream transitions(outdir + "/transitions_" + name + "_" + std::to_string(hours) + ".tsv");
        metrics << "retained\tprimary_retained\tuncertain\tfailure\trho\tquartile_change\n"; std::int64_t transition[4][6]{}; std::uniform_real_distribution<double> phase_dist(0.0, hours*60.0);
        for (int r=0;r<replicates;++r) {
            double phase=phase_dist(rng); int retained=0, primary=0, uncertain=0, failure=0; std::vector<int> hfail(hospitals,0);
            for (const auto& e: episodes) { int state=classify(thin(e.values,phase,hours*60.0,rule),e.end); int transition_state=state<0?0:(state==4?1:state+2); ++transition[e.category][transition_state]; bool bad = state<0 || state==4 || state==2 || state==3; retained += state>=0; primary += state>=0 && state!=4; uncertain += state==2 || state==3; failure += bad; hfail[e.hospital]+=bad; }
            std::vector<double> thinned; for (int h:selected) thinned.push_back(static_cast<double>(hfail[h])/hcount[h]); auto thinrank=ranks(thinned); int changed=0; for (size_t i=0;i<selected.size();++i) { int a=std::min(4,std::max(1,static_cast<int>(std::ceil(4*rawrank[i]/selected.size())))); int b=std::min(4,std::max(1,static_cast<int>(std::ceil(4*thinrank[i]/selected.size())))); changed += a!=b; }
            metrics << static_cast<double>(retained)/episodes.size() << '\t' << static_cast<double>(primary)/episodes.size() << '\t' << (primary?static_cast<double>(uncertain)/primary:0) << '\t' << static_cast<double>(failure)/episodes.size() << '\t' << correlation(rawrank,thinrank) << '\t' << static_cast<double>(changed)/selected.size() << '\n';
        }
        transitions << "reference_category\tthinned_state\tcount\n"; for(int a=0;a<4;++a) for(int b=0;b<6;++b) transitions << a << '\t' << b << '\t' << transition[a][b] << '\n';
    }
    std::ofstream meta(outdir + "/simulation_meta.tsv"); meta << "episodes\t" << episodes.size() << "\nhospitals\t" << hospitals << "\nrank_hospitals\t" << selected.size() << "\n"; for(int c=0;c<4;++c) meta << "reference_category_" << c << '\t' << refcat[c] << '\n';
    return 0;
}
