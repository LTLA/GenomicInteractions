#include "Rcpp.h"

#include "index_map.h"
#include <stdexcept>
#include <deque>

// [[Rcpp::export(rng=false)]]
Rcpp::List expand_1D_hits(Rcpp::IntegerVector query_hits, Rcpp::IntegerVector subject_hits,
    Rcpp::IntegerVector query_indices, Rcpp::IntegerVector query_order,
    Rcpp::IntegerVector subject_indices, Rcpp::IntegerVector subject_order) 
{
    // Building the maps from regions->interactions.
    if (query_indices.size()!=query_order.size() || 
        subject_indices.size()!=subject_order.size()) 
    {
        throw std::runtime_error("index and order vectors should be the same length");
    }
    index_map qreg_2_int, sreg_2_int;
    fill_map(qreg_2_int, query_indices.begin(), query_indices.end(), query_order.begin());
    fill_map(sreg_2_int, subject_indices.begin(), subject_indices.end(), subject_order.begin());

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
            sout.insert(sout.end(), sptr, sptr+sn);
        }
    }

    return Rcpp::List::create(
        Rcpp::IntegerVector(qout.begin(), qout.end()),
        Rcpp::IntegerVector(sout.begin(), sout.end())
    );
}
