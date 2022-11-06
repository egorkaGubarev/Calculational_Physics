#include <algorithm>
#include <iostream>
#include <cmath>
#include <fstream>
#include <list>
#include <vector>

typedef unsigned int uint;
typedef float my_type;

template <typename type>
struct Point
{
    type x_, y_, z_;

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
    std::vector<std::reference_wrapper<Node<type>>> connected_nodes_;

    Node(Point<type> point, const bool is_rigid) : point_{ point }, is_rigid_{ is_rigid } {}

    bool is_connected(Node<type>& other)
    {
        Node<type>* other_p = &other;
        
        for (Node<type>& connected : connected_nodes_) {
            Node<type>* connected_p = &connected;

            if (other_p == connected_p) {
                return true;
            }
        }

        return false;
    }

    void connect(Node<type>& other)
    {
        if (!is_connected(other)) {
            connected_nodes_.push_back(other);
            std::vector<std::reference_wrapper<Node<type>>>& connections = other.connected_nodes_;
            connections.push_back(*this);
        }
    }

    void disconnect(Node<type>& other)
    {
        Node<type>* to_disconnect_p = &other;
        const size_t connections = std::size(connected_nodes_);
        for (uint i = 0; i < connections; ++i) {
            Node<type>& node = connected_nodes_[i];
            Node<type>* node_p = &node;
            if (node_p == to_disconnect_p) {
                std::cout << "Disconnect" << '\n';
            }
        }
    }
};

template <typename type>
struct Segment
{
    type length_;
    Node<type>& begin_, end_;

    Segment(Node<type>& begin, Node<type>& end) : begin_{ begin }, end_{ end }
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

        length_ = (type) std::sqrt(std::pow(delta_x, 2) + std::pow(delta_y, 2) + std::pow(delta_y, 2));

    }
    
    Node<type> get_center() const
    {
        const Point<type>& begin_point = begin_.point_;
        const Point<type>& end_point = end_.point_;

        const type x_begin = begin_point.x_;
        const type y_begin = begin_point.y_;
        const type z_begin = begin_point.z_;

        const type x_end = end_point.x_;
        const type y_end = end_point.y_;
        const type z_end = end_point.z_;

        const type x_center = (x_begin + x_end) / 2;
        const type y_center = (y_begin + y_end) / 2;
        const type z_center = (z_begin + z_end) / 2;

        Point<type> center(x_center, y_center, z_center);
        Node<type> center_node(center, false);
        return center_node;
    }

    ~Segment() = default;
};

template <typename type>
struct Triangle
{
    type length_, square_, max_length_, min_length_;
    Node<type>& node_1_, node_2_, node_3_;

    Triangle(Node<type>& node_1, Node<type>& node_2, Node<type>& node_3):
    node_1_{ node_1}, node_2_{ node_2}, node_3_{ node_3}
    {
        const Segment<type> first(node_1_, node_2_);
        const Segment<type> second(node_1_, node_3_);
        const Segment<type> theird(node_2_, node_3_);

        const type length_1 = first.length_;
        const type length_2 = second.length_;
        const type length_3 = theird.length_;

        max_length_ = std::max(length_1, std::max(length_2, length_3));
        min_length_ = std::min(length_1, std::min(length_2, length_3));

        length_ = length_1 + length_2 + length_3;
        const type half_length = length_ / 2;
        square_ = std::sqrt(half_length * (half_length - length_1) * (half_length - length_2) * (half_length - length_3));
    }

    bool is_equal(const Triangle<type> other) const
    {
        const Node<type>& other_node_1 = other.node_1_;
        const Node<type>& other_node_2 = other.node_2_;
        const Node<type>& other_node_3 = other.node_3_;

        const Node <type>* node_1_ptr = &node_1_;
        const Node <type>* node_2_ptr = &node_2_;
        const Node <type>* node_3_ptr = &node_3_;

        const Node<type>* other_node_1_ptr = &other_node_1;
        const Node<type>* other_node_2_ptr = &other_node_2;
        const Node<type>* other_node_3_ptr = &other_node_3;

        if (node_1_ptr == other_node_1_ptr && node_2_ptr == other_node_2_ptr && node_3_ptr == other_node_3_ptr) {
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
    my_type x_max_, x_min_, y_max_, y_min_, z_max_, z_min_, x_center_, y_center_, z_center_;
    std::list<Node<type>> soft_nodes_, rigid_nodes_;
    std::list<Triangle<type>> triangles_;

    Circuit(Point<type> generator(const type t), const uint dots)
    {
        const type step = (type) 1 / dots;
        Point<type> zero = generator(0);
        create_rigid_node(zero);

        for (uint dot = 1; dot < dots; ++ dot) {
            const type t = step * dot;
            Point<type> new_point = generator(t);
            create_rigid_node(new_point);
            Node<type>& prev = get_last_but_one_rigid_node();
            Node<type>& current = get_last_rigid_node();
            connect(prev, current);
        }

        connect_last_and_zero();
    }

    void connect(Node<type>& first, Node<type>& second)
    {
        first.connect(second);
        find_new_triangles(first, second);
    }

    void connect_last_and_zero()
    {
        Node<type>& last = get_last_rigid_node();
        typename std::list<Node<type>>::iterator begin = std::begin(rigid_nodes_);
        Node<type>& first = *begin;
        connect(first, last);
    }

    // Use when only primary node is created
    void create_primary_surf()
    {
        typename std::list<Node<type>>::iterator begin = std::begin(soft_nodes_);
        Node<type>& center = *begin;

        for (Node<type>& to_connect : rigid_nodes_) {
            connect(center, to_connect);
        }
    }

    void create_primary_node()
    {
        Point<type> center(x_center_, y_center_, z_center_);
        create_soft_node(center);
    }

    void create_rigid_node(const Point<type> point)
    {
        Node<type> node(point, true);
        rigid_nodes_.push_back(node);
    }

    void create_soft_node(const Point<type> point)
    {
        Node<type> node(point, false);
        soft_nodes_.push_back(node);
    }

    void find_new_triangles(Node<type>& first, Node<type>& second)
    {
        const std::vector<std::reference_wrapper<Node<type>>>& first_connect = first.connected_nodes_;
        const std::vector<std::reference_wrapper<Node<type>>>& second_connect = second.connected_nodes_;

        const Node<type>* const first_ptr = &first;
        const Node<type>* const second_ptr = &second;

        for (Node<type>& theird : first_connect) {
            Node<type>* const theird_ptr = &theird;

            if (theird_ptr != second_ptr) {

                for (const Node<type>& node_connect_with_second : second_connect) {
                    const Node<type>* const node_connect_with_second_ptr = &node_connect_with_second;

                    if (theird_ptr == node_connect_with_second_ptr) {
                        Triangle<type> new_triangle(first, second, theird);
                        const bool is_old = is_triangle_exist(new_triangle);
                        if (!is_old) {
                            triangles_.push_back(new_triangle);
                        }
                    }
                }
            }
        }
    }

    Triangle<type>& find_max_length_triangle()
    {
        bool found = false;
        type length = 0;
        Triangle<type>* max_triangle_ptr = nullptr;
        for (Triangle<type>& triangle : triangles_) {
            if (triangle.length_ > length || !found) {
                found = true;
                length = triangle.length_;
                max_triangle_ptr = &triangle;
            }
        }
        return *max_triangle_ptr;
    }

    Triangle<type>& find_min_length_triangle()
    {
        bool found = false;
        type length = 0;
        Triangle<type>* min_triangle_ptr = nullptr;
        for (Triangle<type>& triangle : triangles_) {
            if (triangle.length_ < length || !found) {
                found = true;
                length = triangle.length_;
                min_triangle_ptr = &triangle;
            }
        }
        return *min_triangle_ptr;
    }

    Triangle<type>& find_max_square_triangle()
    {
        bool found = false;
        type square = 0;
        Triangle<type>* max_triangle_ptr = nullptr;
        for (Triangle<type>& triangle : triangles_) {
            if (triangle.square_ > square || !found) {
                found = true;
                square = triangle.square_;
                max_triangle_ptr = &triangle;
            }
        }
        return *max_triangle_ptr;
    }

    Triangle<type>& find_min_square_triangle()
    {
        bool found = false;
        type square = 0;
        Triangle<type>* min_triangle_ptr = nullptr;
        for (Triangle<type>& triangle : triangles_) {
            if (triangle.square_ < square || !found) {
                found = true;
                square = triangle.square_;
                min_triangle_ptr = &triangle;
            }
        }
        return *min_triangle_ptr;
    }

    void get_borders_coords()
    {
        bool x_max_found = false, x_min_found = false, y_max_found = false, y_min_found = false; bool z_min_found = false, z_max_found = false;

        for (const Node<my_type>& node : rigid_nodes_) {
            const Point<my_type>& point = node.point_;

            const my_type x = point.x_;
            const my_type y = point.y_;
            const my_type z = point.z_;

            if (!x_max_found || x > x_max_) {
                x_max_ = x;
                x_max_found = true;
            }

            if (!x_min_found || x < x_min_) {
                x_min_ = x;
                x_min_found = true;
            }

            if (!y_max_found || y > y_max_) {
                y_max_ = y;
                y_max_found = true;
            }

            if (!y_min_found || y < y_min_) {
                y_min_ = y;
                y_min_found = true;
            }

            if (!z_max_found || z > z_max_) {
                z_max_ = z;
                z_max_found = true;
            }

            if (!z_min_found || z < z_min_) {
                z_min_ = z;
                z_min_found = true;
            }
        }
    }

    void get_center()
    {
        x_center_ = (x_min_ + x_max_) / 2;
        y_center_ = (y_min_ + y_max_) / 2;
        z_center_ = (z_min_ + z_max_) / 2;
    }

    // Use when you have at least two nodes only
    Node<type>& get_last_but_one_rigid_node()
    {
        typename std::list<Node<type>>::iterator end = std::end(rigid_nodes_);
        typename std::list<Node<type>>::iterator last = std::prev(end, 2);
        Node<type>& node = *last;
        return node;
    }

    // Use when you have at least one node only
    Node<type>& get_last_rigid_node()
    {
        typename std::list<Node<type>>::iterator end = std::end(rigid_nodes_);
        typename std::list<Node<type>>::iterator last = std::prev(end);
        Node<type>& node = *last;
        return node;
    }

    void get_surf()
    {
        get_borders_coords();
        get_center();
        create_primary_node();
        create_primary_surf();
    }

    bool is_triangle_exist(const Triangle<type>& triangle)
    {
        for (const Triangle<type>& other : triangles_) {
            if (triangle.is_equal(other)) {
                return true;
            }
        }
        return false;
    }

    std::ofstream& log(std::ofstream& out) const
    {
        log_rigid(out);
        log_soft(out);
        log_connect(out);

        return out;
    }

    std::ofstream& log_connect(std::ofstream& out) const
    {
        for (const Node<type>& node : soft_nodes_) {
            const Point<type>& point = node.point_;
            const std::vector<std::reference_wrapper<Node<type>>>& connect = node.connected_nodes_;

            const type x = point.x_;
            const type y = point.y_;
            const type z = point.z_;

            for (const Node<type>& connect_node : connect) {
                const Point<type>& connect_point = connect_node.point_;

                const type connect_x = connect_point.x_;
                const type connect_y = connect_point.y_;
                const type connect_z = connect_point.z_;

                out << x << ' ' << y << ' ' << z << '\n';
                out << connect_x << ' ' << connect_y << ' ' << connect_z << '\n';
            }
        }

        return out;
    }

    std::ofstream& log_rigid(std::ofstream& out) const
    {
        for (const Node<type>& node : rigid_nodes_) {
            const Point<type>& point = node.point_;

            const type x = point.x_;
            const type y = point.y_;
            const type z = point.z_;

            out << x << ' ' << y << ' ' << z << '\n';
        }

        return out;
    }

    std::ofstream& log_soft(std::ofstream& out) const
    {
        for (const Node<type>& node : soft_nodes_) {
            const Point<type>& point = node.point_;

            const type x = point.x_;
            const type y = point.y_;
            const type z = point.z_;

            out << x << ' ' << y << ' ' << z << '\n';
        }

        return out;
    }

    bool partition_done()
    {
        const Triangle<type>& max_square_triangle = find_max_square_triangle();
        const Triangle<type>& min_square_triangle = find_min_square_triangle();
        const Triangle<type>& max_length_triangle = find_max_length_triangle();
        const Triangle<type>& min_length_triangle = find_min_length_triangle();

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

    ~Circuit() = default;
};

template <typename type>
Point<type> circle(const type t)
{
    const type pi = (type) std::acos(-1);
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
    const type vertic_amplit = (type) 1;
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
        y = (type) (1 - std::cos(angle));
    }
    else {
        y = (type) (std::cos(angle) - 1);
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
    std::vector<std::reference_wrapper<std::ofstream>> files{file_circle, file_deform_circle, file_eight, file_cardioid};
    std::vector<std::reference_wrapper<Circuit<my_type>>> circs{ circ_circle, circ_deform_circle, circ_eight, circ_cardioid };
    const size_t circs_amount = std::size(circs);

    for (size_t circ_numb = 0; circ_numb < circs_amount; ++circ_numb) {
        Circuit<my_type>& circ = circs[circ_numb];
        circ.get_surf();
        std::ofstream& file = files[circ_numb];
        circ.log(file);
    }
    
    for (std::ofstream& file : files) {
        file.close();
    }

    const bool done = circ_cardioid.partition_done();

    return 0;
}