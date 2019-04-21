#include <unordered_map>

typedef std::unordered_map<int, std::pair<const int*, int> > index_map;

template<class IT>
void fill_map(index_map& map, IT start, IT end, const int* subject) {
    if (start==end) {
        return;
    }

    int n=1;
    int previous=*start;
    ++start;
    ++subject;

    while (start!=end) {
        if (*start!=previous) {
            if (*start < previous) {
                throw std::runtime_error("indices should be sorted");
            }
            map[previous]=std::make_pair(subject - n, n);
            previous=*start;
            n=1;
        } else {
            ++n;
        }
        ++start;
        ++subject;
    }

    map[previous]=std::make_pair(subject - n, n);
    return;

}
