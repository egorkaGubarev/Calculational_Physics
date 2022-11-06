#include <array>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <unordered_set>

typedef unsigned int uint;
typedef float my_type;

template <typename type>
struct Point
{
    type x_, y_, z_;

    Point(): x_{0}, y_{0}, z_{} {}

    Point(const type x, const type y, const type z) : x_{ x }, y_{ y }, z_{ z } {}

    std::ofstream& Log(std::ofstream& out) const
    {
        out << x_ << ' ' << y_ << ' ' << z_ << '\n';
        return out;
    }

    ~Point() = default;
};

template <typename type>
struct Node
{
    bool is_rigid_;
    Point<type> point_;
    std::unordered_set<uint> connected_nodes_;

    Node(): is_rigid_{true}
    {
        point_ = Point<type>();
    }

    Node(Point<type> point, const bool is_rigid) : point_{ point }, is_rigid_{ is_rigid } {}

    type get_x() const
    {
        return point_.x_;
    }

    type get_y() const
    {
        return point_.y_;
    }

    type get_z() const
    {
        return point_.z_;
    }

    bool is_connected(const uint number)
    {
        if (connected_nodes_.contains(number)) {
            return true;
        }
        else {
            return false;
        }
    }

    void connect(const uint number)
    {
        if (!is_connected(number)) {
            connected_nodes_.insert(number);
        }
    }

    void disconnect(const uint number)
    {
        connected_nodes_.erase(number);
    }
};

template <typename type>
struct Segment
{
    type length_;
    const Node<type>& begin_, end_;

    Segment(const Node<type>& begin, const Node<type>& end) : begin_{ begin }, end_{ end }
    {
        const Point<type>& point_begin = begin_.point_;
        const Point<type>& point_end = end_.point_;

        const type x_begin = point_begin.x_;
        const type y_begin = point_begin.y_;
        const type z_begin = point_begin.z_;

        const type x_end = point_end.x_;
        const type y_end = point_end.y_;
        const type z_end = point_end.z_;

        const type delta_x = x_end - x_begin;
        const type delta_y = y_end - y_begin;
        const type delta_z = z_end - z_begin;

        length_ = (type)std::sqrt(std::pow(delta_x, 2) + std::pow(delta_y, 2) + std::pow(delta_y, 2));

    }

    ~Segment() = default;
};

template <typename type>
struct Triangle
{
    type length_, square_, max_length_, min_length_;
    uint node_1_, node_2_, node_3_;

    Triangle(const uint node_1_number, const uint node_2_number, const uint node_3_number,
        const Node<type>& node_1, const Node<type>& node_2, const Node<type>& node_3)
        : node_1_{ node_1_number }, node_2_{ node_2_number }, node_3_{ node_3_number }
    {
        const Segment<type> first(node_1, node_2);
        const Segment<type> second(node_1, node_3);
        const Segment<type> theird(node_2, node_3);

        const type length_1 = first.length_;
        const type length_2 = second.length_;
        const type length_3 = theird.length_;

        max_length_ = std::max(length_1, std::max(length_2, length_3));
        min_length_ = std::min(length_1, std::min(length_2, length_3));

        length_ = length_1 + length_2 + length_3;
        const type half_length = length_ / 2;
        square_ = std::sqrt(half_length * (half_length - length_1) * (half_length - length_2) * (half_length - length_3));
    }

    bool is_equal(const Triangle<type>& other) const
    {
        const uint other_node_1 = other.node_1_;
        const uint other_node_2 = other.node_2_;
        const uint other_node_3 = other.node_3_;

        if ((node_1_ == other_node_1) && (node_2_ == other_node_2) && (node_3_ == other_node_3)) {
            return true;
        }
        else {
            return false;
        }
    }

    ~Triangle() = default;
};

template <typename type>
struct Circuit
{
    const uint rigid_nodes_;
    std::vector<Node<type>> nodes_;
    std::vector<Triangle<type>> triangles_;

    Circuit(Point<type> generator(const type t), const uint dots) : rigid_nodes_{dots}
    {
        const type step = (type)1 / dots;
        Point<type> zero = generator(0);
        Node<type> zero_node(zero, true);
        nodes_.resize(dots);
        nodes_[0] = zero_node;

        for (uint dot = 1; dot < dots; ++dot) {
            const uint prev_dot_number = dot - 1;
            const type t = step * dot;
            Point<type> new_point = generator(t);
            Node<type> new_node(new_point, true);
            nodes_[dot] = new_node;
            connect(prev_dot_number, dot);
        }

        const uint last_index = dots - 1;
        connect(0, last_index);
    }

    void connect(const uint first, const uint second)
    {
        Node<type>& first_node = nodes_[first];
        Node<type>& second_node = nodes_[second];

        first_node.connect(second);
        second_node.connect(first);

        find_new_triangles(first, second);
    }

    void connect(const uint node, const Triangle<type>& triangle)
    {
        uint node_1 = triangle.node_1_;
        uint node_2 = triangle.node_2_;
        uint node_3 = triangle.node_3_;

        connect(node, node_1);
        connect(node, node_2);
        connect(node, node_3);
    }

    void create_primary_surf(const uint soft_node)
    {
        for (uint node = 0; node < rigid_nodes_; ++ node) {
            connect(soft_node, node);
        }
    }

    void create_primary_node()
    {
        const std::array<type, 3> mean = get_mean_coords();
        Point<type> center(mean[0], mean[1], mean[2]);
        Node<type> center_node(center, false);
        nodes_.push_back(center_node);
    }

    void erase_triangle(const uint number)
    {
        typename std::vector<Triangle<type>>::iterator to_erase = std::begin(triangles_);
        std::advance(to_erase, number);
        triangles_.erase(to_erase);
    }

    void find_new_triangles(const uint first, const uint second)
    {
        const Node<type>& first_node = nodes_[first];
        const Node<type>& second_node = nodes_[second];

        const std::unordered_set<uint>& first_connect = first_node.connected_nodes_;
        const std::unordered_set<uint>& second_connect = second_node.connected_nodes_;

        for (const uint theird : first_connect) {
            if (theird != second) {
                for (const uint node_connect_with_second : second_connect) {
                    if (theird == node_connect_with_second) {
                        const Node<type>& third_node = nodes_[theird];
                        Triangle<type> new_triangle(first, second, theird, first_node, second_node, third_node);

                        if (!(triangle_exist(new_triangle))) {
                            triangles_.push_back(new_triangle);
                        }

                    }
                }
            }
        }
    }

    uint find_max_length_triangle_number() const
    {
        bool found = false;
        uint max_triangle_number = 0;
        type length = 0;
        const size_t triangles = std::size(triangles_);
        for (uint triangle_number = 0; triangle_number < triangles; ++ triangle_number) {
            const Triangle<type>& triangle = triangles_[triangle_number];
            if (triangle.length_ > length || !found) {
                found = true;
                length = triangle.length_;
                max_triangle_number = triangle_number;
            }
        }
        return max_triangle_number;
    }

    uint find_min_length_triangle_number() const
    {
        bool found = false;
        uint min_triangle_number = 0;
        type length = 0;
        const size_t triangles = std::size(triangles_);
        for (uint triangle_number = 0; triangle_number < triangles; ++triangle_number) {
            const Triangle<type>& triangle = triangles_[triangle_number];
            if (triangle.length_ < length || !found) {
                found = true;
                length = triangle.length_;
                min_triangle_number = triangle_number;
            }
        }
        return min_triangle_number;
    }

    uint find_max_square_triangle_number() const
    {
        bool found = false;
        uint max_triangle_number = 0;
        type square = 0;
        const size_t triangles = std::size(triangles_);
        for (uint triangle_number = 0; triangle_number < triangles; ++triangle_number) {
            const Triangle<type>& triangle = triangles_[triangle_number];
            if (triangle.square_ > square || !found) {
                found = true;
                square = triangle.length_;
                max_triangle_number = triangle_number;
            }
        }
        return max_triangle_number;
    }

    uint find_min_square_triangle_number() const
    {
        bool found = false;
        uint min_triangle_number = 0;
        type square = 0;
        const size_t triangles = std::size(triangles_);
        for (uint triangle_number = 0; triangle_number < triangles; ++triangle_number) {
            const Triangle<type>& triangle = triangles_[triangle_number];
            if (triangle.square_ < square || !found) {
                found = true;
                square = triangle.length_;
                min_triangle_number = triangle_number;
            }
        }
        return min_triangle_number;
    }

    void improve_partition()
    {
        if (true) {
            const uint max_triangle_number = find_max_square_triangle_number();
            const Triangle<type>& triangle = triangles_[max_triangle_number];

            const uint node_number_1 = triangle.node_1_;
            const uint node_number_2 = triangle.node_2_;
            const uint node_number_3 = triangle.node_3_;

            const Node<type>& node_1 = nodes_[node_number_1];
            const Node<type>& node_2 = nodes_[node_number_2];
            const Node<type>& node_3 = nodes_[node_number_3];

            const Point<type>& center = get_mass_center(node_1, node_2, node_3);
            Node<type> node(center, true);
            nodes_.push_back(node);
            const uint node_number = uint(std::size(nodes_) - 1);
            connect(node_number, triangle);
            erase_triangle(max_triangle_number);
        }
    }

    std::array<type, 3> get_mean_coords() const
    {
        type x_sum = 0;
        type y_sum = 0;
        type z_sum = 0;

        for (uint node_number = 0; node_number < rigid_nodes_; ++node_number) {
            const Node<type>& node = nodes_[node_number];

            const type x = node.get_x();
            const type y = node.get_y();
            const type z = node.get_z();

            x_sum += x;
            y_sum += y;
            z_sum += z;
        }

        const type x_mean = x_sum / rigid_nodes_;
        const type y_mean = y_sum / rigid_nodes_;
        const type z_mean = z_sum / rigid_nodes_;

        const std::array<type, 3> mean{ x_mean, y_mean, z_mean };
        return mean;
    }

    void get_surf()
    {
        create_primary_node();
        create_primary_surf(rigid_nodes_);
    }

    std::ofstream& log(std::ofstream& out) const
    {
        const size_t nodes = std::size(nodes_);
        const size_t soft_nodes = nodes - rigid_nodes_;
        out << rigid_nodes_ << ' ' << soft_nodes << ' ' << 0 << '\n';

        log_rigid(out);
        log_soft(out);
        log_connect(out);

        return out;
    }

    std::ofstream& log_connect(std::ofstream& out) const
    {
        const size_t nodes = std::size(nodes_);
        for (size_t node_number = rigid_nodes_; node_number < nodes; ++ node_number) {
            const Node<type> node = nodes_[node_number];
            const std::unordered_set<uint>& connect = node.connected_nodes_;

            const type x = node.get_x();
            const type y = node.get_y();
            const type z = node.get_z();

            for (uint connect_node_number : connect) {
                const Node<type> connect_node = nodes_[connect_node_number];

                const type connect_x = connect_node.get_x();
                const type connect_y = connect_node.get_y();
                const type connect_z = connect_node.get_z();

                out << x << ' ' << y << ' ' << z << '\n';
                out << connect_x << ' ' << connect_y << ' ' << connect_z << '\n';
            }
        }

        return out;
    }

    std::ofstream& log_rigid(std::ofstream& out) const
    {
        for (uint node_number = 0; node_number < rigid_nodes_; ++ node_number) {
            const Node<type>& node = nodes_[node_number];

            const type x = node.get_x();
            const type y = node.get_y();
            const type z = node.get_z();

            out << x << ' ' << y << ' ' << z << '\n';
        }

        return out;
    }

    std::ofstream& log_soft(std::ofstream& out) const
    {
        const size_t nodes = std::size(nodes_);
        for (uint node_number = rigid_nodes_; node_number < nodes; ++ node_number) {
            const Node<type>& node = nodes_[node_number];

            const type x = node.get_x();
            const type y = node.get_y();
            const type z = node.get_z();

            out << x << ' ' << y << ' ' << z << '\n';
        }

        return out;
    }

    bool partition_done()
    {
        const uint max_square_triangle_number = find_max_square_triangle_number();
        const uint min_square_triangle_number = find_min_square_triangle_number();
        const uint max_length_triangle_number = find_max_length_triangle_number();
        const uint min_length_triangle_number = find_min_length_triangle_number();

        const Triangle<type> max_square_triangle = triangles_[max_square_triangle_number];
        const Triangle<type> min_square_triangle = triangles_[min_square_triangle_number];
        const Triangle<type> max_length_triangle = triangles_[max_length_triangle_number];
        const Triangle<type> min_length_triangle = triangles_[min_length_triangle_number];

        const type max_square = max_square_triangle.square_;
        const type min_square = min_square_triangle.square_;
        const type max_lenght = max_length_triangle.length_;
        const type min_length = min_length_triangle.length_;

        if (max_square / min_square > 3 || max_lenght / min_length > std::sqrt(3)) {
            return false;
        }
        else {
            return true;
        }
    }

    bool triangle_exist(const Triangle<type>& triangle) {
        for (const Triangle<type>& other : triangles_) {
            if (triangle.is_equal(other)) {
                return true;
            }
        }
        return false;
    }

    ~Circuit() = default;
};

template <typename type>
Point<type> get_mass_center(const Node<type>& node_1, const Node<type>& node_2, const Node<type>& node_3)
{
    const type x_1 = node_1.get_x();
    const type x_2 = node_2.get_x();
    const type x_3 = node_3.get_x();

    const type y_1 = node_1.get_y();
    const type y_2 = node_2.get_y();
    const type y_3 = node_3.get_y();

    const type z_1 = node_1.get_z();
    const type z_2 = node_2.get_z();
    const type z_3 = node_3.get_z();

    const type x = (x_1 + x_2 + x_3) / 3;
    const type y = (y_1 + y_2 + y_3) / 3;
    const type z = (z_1 + z_2 + z_3) / 3;

    Point<type> point(x, y, z);
    return point;
}

template <typename type>
Point<type> circle(const type t)
{
    const type pi = (type)std::acos(-1);
    const type angle = 2 * pi * t;

    const type x = std::cos(angle);
    const type y = std::sin(angle);
    const type z = 0;

    const Point<type> point(x, y, z);
    return point;
}

template <typename type>
Point<type> deform_circle(const type t)
{
    const type vertic_amplit = (type)1;
    const type vertic_frequency = 10;

    const type pi = (type)std::acos(-1);
    const type angle = 2 * pi * t;

    const type x = std::cos(angle);
    const type y = std::sin(angle);
    const type z = vertic_amplit * std::sin(angle * vertic_frequency);

    const Point<type> point(x, y, z);
    return point;
}

template <typename type>
Point<type> eight(const type t)
{
    const type vertic_amplit = (type)0.1;
    const type vertic_frequency = 0.25;

    const type pi = (type)std::acos(-1);
    const type angle = 4 * pi * t;
    type y = 0;
    const type x = (type)std::sin(angle);

    if (angle < 2 * pi) {
        y = (type)(1 - std::cos(angle));
    }
    else {
        y = (type)(std::cos(angle) - 1);
    }

    const type z = vertic_amplit * std::sin(angle * vertic_frequency);
    const Point<type> point(x, y, z);
    return point;
}

template <typename type>
Point<type> cardioid(const type t)
{
    const type vertic_amplit = 1;
    const type vertic_frequency = 1;

    const type pi = (type)std::acos(-1);
    const type angle = 2 * pi * t;

    const type x = (type)(2 * std::cos(angle) - std::cos(2 * angle));
    const type y = (type)(2 * std::sin(angle) - std::sin(2 * angle));
    const type z = vertic_amplit * std::sin(angle * vertic_frequency);

    const Point<type> point(x, y, z);
    return point;
}

int main()
{
    const uint dots_on_circ = 10;
    std::ofstream aux, file_circle, file_deform_circle, file_eight, file_cardioid;

    aux.open("aux_info.txt");
    file_circle.open("circle.txt");
    file_deform_circle.open("deform_circle.txt");
    file_eight.open("eight.txt");
    file_cardioid.open("cardioid.txt");

    Circuit<my_type> circ_circle(circle, dots_on_circ), circ_deform_circle(deform_circle, dots_on_circ), circ_eight(eight, dots_on_circ), circ_cardioid(cardioid, dots_on_circ);

    aux << dots_on_circ << '\n';
    aux.close();
    std::vector<std::reference_wrapper<std::ofstream>> files{ file_circle, file_deform_circle, file_eight, file_cardioid };
    std::vector<std::reference_wrapper<Circuit<my_type>>> circs{ circ_circle, circ_deform_circle, circ_eight, circ_cardioid };
    const size_t circs_amount = std::size(circs);

    for (size_t circ_numb = 0; circ_numb < circs_amount; ++circ_numb) {
        Circuit<my_type>& circ = circs[circ_numb];
        circ.get_surf();
        for (uint i = 0; i < 10; ++i) {
            circ.improve_partition();
        }
        std::ofstream& file = files[circ_numb];
        circ.log(file);
    }

    for (std::ofstream& file : files) {
        file.close();
    }

    return 0;
}
