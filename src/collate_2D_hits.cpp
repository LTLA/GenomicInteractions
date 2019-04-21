#include "Rcpp.h"

#include "index_map.h"
#include <unordered_map>
#include <stdexcept>
#include <deque>

/* This function is rather difficult to grasp as there are multiple levels of indirection:
 *   1. Query anchor indices to query regions (specifying the interactions in GenomicInteractions)
 *   2. Query regions to subject regions (from findOverlaps between genomic regions)
 *   3. Subject regions to subject anchor indices (going backwards to get to the interactions)
 *   4. Subject anchor indices to subject interaction index (to get the actual index of the interactions)
 *
 * To use the left anchor region as an example:
 *   - 'query_indices_left' is of length equal to the number of query interactions,
 *     and contains the left anchor index for each query interaction.
 *     This is relevant to Indirection #1.
 *   - 'query_hits_left' and 'subject_hits_left' are the queryHits() and subjectHits() of the Hits object 
 *     that results from a findOverlaps() on the left query regions and left subject regions.
 *     This is relevant to Indirection #2.
 *   - 'subject_indices_left' is of length equal to the number of subject interactions,
 *     and contains the left anchor index for each subject interaction.
 *     This is relevant to Indirection #3.
 *   - 'subject_order_left' is of length equal to the number of subject interactions,
 *     and contains the permutation index for each subject interaction to ensure that 'subject_indices_left' is sorted.
 *     This is relevant to Indirection #4.
 * 
 * The same applies for the right anchor regions, 
 * which can be effectively considered in isolation from the left for the most part.
 * (The two only come together at the end when intersecting interaction indices.)
 */

// [[Rcpp::export(rng=false)]]
Rcpp::List collate_2D_hits(Rcpp::IntegerVector query_indices_left, Rcpp::IntegerVector query_indices_right,
    Rcpp::IntegerVector query_hits_left, Rcpp::IntegerVector subject_hits_left,
    Rcpp::IntegerVector query_hits_right, Rcpp::IntegerVector subject_hits_right,
    Rcpp::IntegerVector subject_indices_left, Rcpp::IntegerVector subject_order_left,
    Rcpp::IntegerVector subject_indices_right, Rcpp::IntegerVector subject_order_right) 
{
    // Creating maps to go from query regions->subject regions.
    if (query_hits_left.size()!=subject_hits_left.size() ||
           query_hits_right.size()!=subject_hits_right.size()) 
    {
        throw std::runtime_error("query and subject hits vectors should be the same length");
    }
    index_map left_reg, right_reg;
    fill_map(left_reg, query_hits_left.begin(), query_hits_left.end(), subject_hits_left.begin());
    fill_map(right_reg, query_hits_right.begin(), query_hits_right.end(), subject_hits_right.begin());

    // Creating maps to go from subject regions -> indices.
    if (subject_indices_left.size()!=subject_order_left.size() ||
           subject_indices_right.size()!=subject_order_right.size()) 
    {
        throw std::runtime_error("subject index and order vectors should be the same length");
    }
    index_map left_ind, right_ind;
    fill_map(left_ind, subject_indices_left.begin(), subject_indices_left.end(), subject_order_left.begin());
    fill_map(right_ind, subject_indices_right.begin(), subject_indices_right.end(), subject_order_right.begin());
    
    // Iterating across all queries.
    if (query_indices_left.size()!=query_indices_right.size()) {
        throw std::runtime_error("query index vectors should be the same length");
    }
    auto qlIt=query_indices_left.begin(), qrIt=query_indices_right.begin();
    const size_t N=query_indices_left.size();

    std::deque<int> qout, sout;
    std::unordered_map<int, bool> collected;

    for (size_t i=0; i<N; ++i, ++qlIt, ++qrIt) {
        auto lregIt=left_reg.find(*qlIt);
        if (lregIt==left_reg.end()) {
            continue;
        }
        auto rregIt=right_reg.find(*qrIt);
        if (rregIt==right_reg.end()) {
            continue;
        }

        // We identify all subject regions that overlap with the left region of the current query.
        const auto& lopt=lregIt->second;
        const int* lptr=lopt.first;
        const int ln=lopt.second;

        for (int j=0; j<ln; ++j, ++lptr) {
            // Second search identifies all subject interactions with
            // a left region among the overlapped subject region.
            auto idxIt=left_ind.find(*lptr);
            if (idxIt!=left_ind.end()) {
                const auto& opt=idxIt->second;
                const int* ptr=opt.first;
                const int n=opt.second;
                for (int k=0; k<n; ++k, ++ptr) {
                    collected[*ptr]=false;
                }
            }
        }

        if (collected.empty()) {
            continue;
        }

        // We next identify all subject regions that overlap with the right region of the current query.
        const auto& ropt=rregIt->second;
        const int* rptr=ropt.first;
        const int rn=ropt.second;

        for (int j=0; j<rn; ++j, ++rptr) {
            // Second search identifies all subject interactions with
            // a right region among the overlapped subject region. We then
            // intersect this with those we found for the left overlaps.
            auto idxIt=right_ind.find(*rptr);
            if (idxIt==right_ind.end()) { continue; }

            const auto& opt=idxIt->second;
            const int* ptr=opt.first;
            const int n=opt.second;
            for (int k=0; k<n; ++k, ++ptr) {
                auto it=collected.find(*ptr);
                if (it!=collected.end() && !(it->second)) {
                    it->second=true; // mark this subject as having been recorded already.
                    qout.push_back(i+1); // get back to 1-based indexing.
                    sout.push_back(*ptr);
                }
            }
        }

        collected.clear();
    }

    return Rcpp::List::create(
        Rcpp::IntegerVector(qout.begin(), qout.end()),
        Rcpp::IntegerVector(sout.begin(), sout.end())
    );
}
