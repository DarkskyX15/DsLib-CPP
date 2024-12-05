
#include <iostream>
#include <string>
#include <stack>
#include "../../Graph.hpp"

using dsl::graph::defines::UpdateStrategy;
using namespace dsl::graph;

struct Person {
    std::string id_card;
    std::string name;
    size_t tel;

    bool operator==(const Person& p) const { return id_card == p.id_card; }
    bool operator!=(const Person& p) const { return !operator==(p); }
};

template<>
struct utils::key_selector<Person> {
    typedef std::string key_type;
    static constexpr const std::string& key(const Person& p) {
        return p.name;
    }
    static constexpr std::string& key(Person& p) {
        return p.name;
    }
    static constexpr std::string&& arg(
        std::string&& ,
        std::string&& s,
        size_t&&
    ) { return static_cast<std::string&&>(s); }
};

int main() {
    SimpleGraph<
        Person, bool, false, size_t,
        HashListStorage<size_t, bool, false>,
        DefaultIndexProvider<Person, size_t>
    > g;

    std::string card1 = "idCardNo1";
    std::string name1 = "FS";

    // emplace node
    g.emplaceNode(card1, name1, 515);
    g.emplaceNode("idCardNo2", "FFFF", 1145);
    g.emplaceNode("idCardNo3", "SA", 515);
    g.emplaceNode("idCardNo4", "DS", 1115);
    g.emplaceNode("idCardNo5", "TH", 5520);

    // add edge
    g.addEdgeByKey(name1, "FFFF")
     .addEdgeByKey("FFFF", "DS")
     .addEdgeByKey("FS", "TH")
     .addEdgeByKey("SA", "TH");

    // walk
    algorithms::BFS(g.const_access(),
        [](const Person& x) -> bool {
            std::cout << "{ ID Card: " << x.id_card
                << " Name: " << x.name
                << " TEL: " << x.tel
                << " }\n";
            return true;
        }
    );

    auto start = g.const_access(g.find("DS"));
    auto dest = g.const_access(g.find("SA"));

    typedef decltype(start) __accessor;
    std::stack<const Person*> buf;
    std::unordered_set<size_t> visited;
    visited.insert(start.raw());

    std::function<bool(__accessor)> _rec = 
    [&buf, &_rec, &dest, &visited](__accessor acc) -> bool {
        if (acc == dest) {
            buf.push(&*acc);
            return true;
        }
        bool needed = false;
        acc.updateAdjacent(UpdateStrategy::forth);
        for (auto [idx, vp, wp]: acc.listForth()) {
            if (!visited.contains(idx)) {
                visited.insert(idx);
                needed |= _rec(acc.next(idx));
                if (needed) {
                    buf.push(&*acc);
                    return true;
                }
            }
        }
        return false;
    };

    auto solved = _rec(start);

    if (!solved) {
        std::cout << "Can not reach!\n" ;
        return 0;
    }

    std::cout << "START";
    while (!buf.empty()) {
        std::cout << " -> " << buf.top()->name;
        buf.pop();
    }

    return 0;
}