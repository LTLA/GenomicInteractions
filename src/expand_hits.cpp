#include "Rcpp.h"

#include <unordered_map>
#include <stdexcept>
#include <deque>

typedef std::unordered_map<int, std::pair<const int*, int> > index_map;

void expand_map(index_map& map, Rcpp::IntegerVector indices) {
    if (indices.size()==0) {
        return;
    }

    int n=1, previous=indices[0];
    const int* it=indices.begin();
    const int* end=it+indices.size();

    while ((++it)!=end) {
        if (*it!=previous) {
            if (*it < previous) {
                throw std::runtime_error("indices should be sorted");
            }
            map[previous]=std::make_pair(it - n, n);
            previous=*it;
            n=1;
        } else {
            ++n;
        }
    }

    map[previous]=std::make_pair(it - n, n);
    return;
}

// [[Rcpp::export(rng=false)]]
Rcpp::List expand_hits(Rcpp::IntegerVector query_hits, Rcpp::IntegerVector subject_hits,
    Rcpp::IntegerVector query_indices, Rcpp::IntegerVector subject_indices) 
{
    // Building the maps from regions->interactions.
    index_map qreg_2_int, sreg_2_int;
    expand_map(qreg_2_int, query_indices);
    expand_map(sreg_2_int, subject_indices);

    // Running through the hits and spawning all combinations.
    if (query_hits.size()!=subject_hits.size()) {
        throw std::runtime_error("'query_hits' and 'subject_hits' should be the same length");
    }

    auto qIt=query_hits.begin(), sIt=subject_hits.begin();
    std::deque<int> qout, sout;

    const size_t N=query_hits.size();
    for (size_t i=0; i<N; ++i, ++qIt, ++sIt) {
        auto qmapIt=qreg_2_int.find(*qIt);
        if (qmapIt==qreg_2_int.end()) {
            continue;
        }
        auto smapIt=sreg_2_int.find(*sIt);
        if (smapIt==sreg_2_int.end()) {
            continue;
        }

        const auto& qopt=qmapIt->second;
        const int* qptr=qopt.first;
        const int qn=qopt.second;

        const auto& sopt=smapIt->second;
        const int* sptr=sopt.first;
        const int sn=sopt.second;

        for (int jq=0; jq<qn; ++jq, ++qptr) {
            qout.insert(qout.end(), sn, *qptr);
            sout.insert(sout.begin(), sptr, sptr+sn);
        }
    }

    return Rcpp::List::create(
        Rcpp::IntegerVector(qout.begin(), qout.end()),
        Rcpp::IntegerVector(sout.begin(), sout.end())
    );
}
