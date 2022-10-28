#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <vector>

typedef unsigned int uint;

template <typename type>
class UniformRandomGenerator
{
public:
    explicit UniformRandomGenerator(long long seed)
    {
        if (seed == -1) {
            seed = time(nullptr);
        }
        engine_ = new std::mt19937_64(seed);
    }

    std::vector<type> create_sequence(const unsigned int n, const type min, const type max)
    {
        std::vector<type> sequence(n);
        for (unsigned int i = 0; i < n; ++i) {
            sequence[i] = get_number(min, max);
        }
        return sequence;
    }

    type get_number(type min, type max)
    {
        std::uniform_int_distribution<type> distr(min, max);
        const type number = distr(*engine_);
        return number;
    }

private:
    std::mt19937_64* engine_;
};

struct Spot
{
    bool free;
    bool is_stop;
    uint pass;

    Spot() : free{ true }, is_stop{ false }, pass{ 0 } {}
};

class Road
{
public:
    const uint length_;
    std::vector<Spot> lane_;
    const uint name_;
    Road* beg_conn;
    Road* end_conn;

    Road(const uint length, const uint name) : length_{ length }, name_{ name }, beg_conn{ nullptr }, end_conn{ nullptr }
    {
        lane_.resize(length_);
        Spot default_spot;
        std::fill(std::begin(lane_), std::end(lane_), default_spot);
    }

    ~Road() = default;

    void AddStop(const uint pos)
    {
        Spot& spot_with_stop = lane_[pos];
        spot_with_stop.is_stop = true;
        spot_with_stop.pass = 0;
    }

    void connect(Road& other)
    {
        end_conn = &other;
        other.beg_conn = this;
    }

    std::ofstream& Log(std::ofstream& file) const
    {
        for (uint i = 0; i < length_; ++i) {
            const Spot& spot = lane_[i];
            const bool taken = !spot.free;
            if (taken) {

                file << name_ * length_ + i << '\n';
            }
        }
        return file;
    }

    std::vector<uint> spots_to_vector(std::vector<uint>& spots) const
    {
        for (uint i = 0; i < length_; ++i) {
            const Spot& spot = lane_[i];
            const bool taken = !spot.free;
            if (taken) {
                uint coordin = name_ * length_ + i;
                spots.push_back(coordin);
            }
        }
        return spots;
    }
};

class Trolley
{
public:
    Trolley(const uint pos, const uint speed, const std::string name, Road& locat) :
        pass_{ 0 }, pass_to_drop_{ 0 }, pos_{ pos }, speed_{ speed }, stat_{ "move" }, name_{ name }, locat_{ &locat }
    {
        TakeSpot();
    }

    ~Trolley() = default;

    void FreeSpot()
    {
        std::vector<Spot>& lane = locat_->lane_;
        Spot& spot = lane[pos_];
        spot.free = true;
    }

    void TakeSpot()
    {
        std::vector<Spot>& lane = locat_->lane_;
        Spot& spot = lane[pos_];
        spot.free = false;
    }

    static bool NeedMoveMore(const uint max_dist, const uint av_dist, const uint dist)
    {
        if (dist < max_dist && dist == av_dist) {
            return true;
        }
        else {
            return false;
        }
    }

    // Use when max_dist <= free_space only
    uint GetRange(const uint max_dist, const Road* locat) const
    {
        uint dist = 1;
        const Road& street = *locat;
        const std::vector<Spot>& lane = street.lane_;
        bool is_free = true;
        while (dist <= max_dist && is_free) {
            const uint pos_to_check = pos_ + dist;
            const Spot& spot_to_check = lane[pos_to_check];
            is_free = spot_to_check.free;
            ++dist;
        }
        if (!is_free) {
            const uint possible_dist = dist - 2;
            return possible_dist;
        }
        else {
            const uint possible_dist = dist - 1;
            return possible_dist;
        }
    }

    uint GetStopRange(const uint max_dist, const Road* locat) const
    {
        uint dist = 1;
        const Road& street = *locat;
        const std::vector<Spot>& lane = street.lane_;
        bool is_stop = false;
        while (dist <= max_dist && !is_stop) {
            const uint pos_to_check = pos_ + dist;
            const Spot& spot_to_check = lane[pos_to_check];
            is_stop = spot_to_check.is_stop;
            ++dist;
        }
        const uint possible_dist = dist - 1;
        return possible_dist;
    }

    // Trolley can't skip street
    void Move(const uint dt)
    {
        if (stat_ == "move") {
            FreeSpot();
            const uint max_dist = speed_ * dt;
            Road& street = *locat_;
            const uint length = street.length_;
            const uint free_space = length - pos_ - 1;
            const uint av_dist = std::min(max_dist, free_space);
            const uint dist = GetRange(av_dist, locat_);
            const uint dist_to_stop = GetStopRange(av_dist, locat_);
            const uint dist_to_move = std::min(dist, dist_to_stop);
            pos_ += dist_to_move;
            if (NeedMoveMore(max_dist, av_dist, dist_to_move)) {
                Road* next_street_p = street.end_conn;
                Road& next_street = *next_street_p;
                std::vector<Spot>& next_lane = next_street.lane_;
                Spot& zero_spot = next_lane[0];
                if (zero_spot.free) {
                    locat_ = next_street_p;
                    pos_ = 0;
                    if (!zero_spot.is_stop) {
                        const uint dist_passed = dist_to_move + 1;
                        if (dist_passed < max_dist) {
                            const uint max_dist_left = max_dist - dist_passed;
                            const uint dist = GetRange(max_dist_left, locat_);
                            const uint dist_to_stop = GetStopRange(av_dist, locat_);
                            const uint dist_to_move = std::min(dist, dist_to_move);
                            pos_ += dist_to_move;
                        }
                    }
                }
            }
            TakeSpot();
        }
    }

    void GetPass(const uint pass_per_sec_board)
    {
        Road& street = *locat_;
        std::vector<Spot>& lane = street.lane_;
        Spot& spot = lane[pos_];
        if (spot.is_stop) {
            if (stat_ == "board") {
                if (spot.pass > 0) {
                    if (spot.pass < pass_per_sec_board) {
                        pass_ += spot.pass;
                        spot.pass = 0;
                        stat_ = "move";
                    }
                    else {
                        spot.pass -= pass_per_sec_board;
                        pass_ += pass_per_sec_board;

                    }
                }
                else {
                    stat_ = "move";
                }
            }
        }
    }

    void DropPass(const uint pass_per_sec_board, UniformRandomGenerator<uint>& generator)
    {
        Road& street = *locat_;
        std::vector<Spot>& lane = street.lane_;
        Spot& spot = lane[pos_];
        if (spot.is_stop) {
            if (stat_ == "move") {
                if (pass_ > 0) {
                    stat_ = "unb";
                    pass_to_drop_ = generator.get_number(0, pass_);
                }
                else {
                    stat_ = "board";
                }
            }
            if (stat_ == "unb") {
                if (pass_to_drop_ > pass_per_sec_board) {
                    pass_ -= pass_per_sec_board;
                    pass_to_drop_ -= pass_per_sec_board;
                }
                else {
                    pass_ -= pass_to_drop_;
                    stat_ = "board";
                }
            }
        }
    }

    void Log() const
    {
        const Road& street = *locat_;
        const uint street_name = street.name_;
        std::cout << name_ << ' ' << street_name << ' ' << pos_ << ' ' << pass_ << '\n';
    }
private:
    uint pass_, pass_to_drop_, pos_, speed_;
    std::string stat_;
    const std::string name_;
    Road* locat_;
};

void push_distance_to_vector(const std::vector<uint>& spots, const uint trolleys, const size_t length, std::vector<float>& vector_dist)
{
    uint prev_corrdin = spots[0];

    for (uint i = 1; i < trolleys; ++i) {
        uint coordin = spots[i];
        const float distance = (float)(coordin - prev_corrdin);
        vector_dist.push_back(distance);
        prev_corrdin = coordin;
    }

    const float circle_distance = (float)(length - spots[trolleys - 1] + spots[0]);
    vector_dist.push_back(circle_distance);
}

float get_entropy(std::vector<float>& vector_dist)
{
    float part_sum = 0;

    for (const float dist : vector_dist) {
        part_sum += (float)(dist * std::log(dist));
    }

    return part_sum;
}

uint get_collapse_time(const float entropy_limit, const uint pass_per)
{
    const uint dt = 1;
    const uint pass_per_sec_board = 1;
    const uint troll = 100;
    const uint length = 2000;
    const uint speed = 10;
    UniformRandomGenerator<uint> generator(-1);

    Road north(length, 0), east(length, 1), south(length, 2), west(length, 3);

    north.AddStop(0);
    east.AddStop(0);
    south.AddStop(0);
    west.AddStop(0);
    north.AddStop(1'000);
    east.AddStop(1'000);
    south.AddStop(1'000);
    west.AddStop(1'000);

    std::vector<std::reference_wrapper<Road>> roads{ north, east, south, west };
    const size_t roads_amount = std::size(roads);
    std::vector<Trolley> trolleys;

    for (uint troll_num = 0; troll_num < troll; ++troll_num) {
        const uint road_num = generator.get_number(0, 3);
        const uint pos = generator.get_number(0, length - 1);
        Trolley new_trolley(pos, speed, "-", roads[road_num]);
        trolleys.push_back(new_trolley);
    }

    std::vector<std::reference_wrapper<Spot>> stops;

    for (Road& street : roads) {
        std::vector<Spot>& lane = street.lane_;
        for (Spot& spot : lane) {
            if (spot.is_stop) {
                stops.push_back(spot);
            }
        }
    }

    north.connect(east);
    east.connect(south);
    south.connect(west);
    west.connect(north);

    uint i = 0;
    float entropy = 2 * entropy_limit;

    while (entropy < entropy_limit) {
        for (Trolley& trolley : trolleys) {
            trolley.DropPass(pass_per_sec_board, generator);
            trolley.GetPass(pass_per_sec_board);
            trolley.Move(dt);
        }

        uint mod = i % pass_per;

        if (mod == 0) {
            for (Spot& stop : stops) {
                const uint new_pass = generator.get_number(0, 1);
                stop.pass += new_pass;
            }
        }

        std::vector<uint> spots;

        for (const Road& street : roads) {
            street.spots_to_vector(spots);
        }

            
        const size_t full_length = length * roads_amount;
        std::vector<float> vector_dist;
        push_distance_to_vector(spots, troll, full_length, vector_dist);

        for (float& dist : vector_dist) {
            dist /= full_length;
        }

        entropy = get_entropy(vector_dist);
        ++i;
    }
        
    return i;
}

int main()
{
    std::ofstream collapse;
    collapse.open("collapse.txt");

    for (uint pass_per = 100; pass_per >= 10; -- pass_per) {
        const uint t = get_collapse_time((float)(-1.2), pass_per);
        collapse << pass_per << ' ' << t << '\n';
        std::cout << "Pass per: " << pass_per << '\n';
        std::cout << "Time: " << t << '\n';
    }
    
    collapse.close();
    return 0;
}
