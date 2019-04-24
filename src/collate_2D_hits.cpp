#include "Rcpp.h"

#include "index_map.h"
#include <unordered_map>
#include <stdexcept>
#include <deque>

/* A convenience class to organize the query region -> subject region -> subject interaction
 * indirection for a given anchor region. This avoids the need to duplicate code for the
 * left and right anchors, as we can simply construct the index_organizer twice.
 */ 

class index_organizer {
public:
    index_organizer(Rcpp::IntegerVector query_hits, Rcpp::IntegerVector subject_hits,
        Rcpp::IntegerVector subject_indices, Rcpp::IntegerVector subject_order) 
    {
        // Going from query region hits to subject region hits.
        if (query_hits.size()!=subject_hits.size()) {
            throw std::runtime_error("query and subject hits vectors should be the same length");
        }
        fill_map(queryReg_to_subjectReg, query_hits.begin(), query_hits.end(), subject_hits.begin());

        // Creating maps to go from subject regions to subject interaction indices.
        if (subject_indices.size()!=subject_order.size()) {
            throw std::runtime_error("subject index and order vectors should be the same length");
        }
        fill_map(subjectReg_to_subjectInt, subject_indices.begin(), subject_indices.end(), subject_order.begin());
        return;
    }

    index_map::const_iterator search_query_region(int i) const {
        return queryReg_to_subjectReg.find(i);
    }

    bool has_query_region(const index_map::const_iterator& it) const {
        return it!=queryReg_to_subjectReg.end();
    }

    index_map::const_iterator search_subject_region(int i) const {
        return subjectReg_to_subjectInt.find(i);
    }

    bool has_subject_region(const index_map::const_iterator& it) const {
        return it!=subjectReg_to_subjectInt.end();
    }
private:
    index_map queryReg_to_subjectReg;
    index_map subjectReg_to_subjectInt;
};

/* A convenience class to organize the results (interaction indices),
 * as well as the intermediate map used to intersect left/right indices.
 * It also helps to avoid nested loops that make it difficult to break 
 * out when skip=true (to avoid a further search for select='arbitrary').
 */ 

class hits_collector {
public:
    hits_collector(bool skip_on_first) : skip(skip_on_first) {}

    void reset_left_store() {
        encountered.clear();
        return;
    };

    bool no_left_hits() const {
        return encountered.empty();
    }

    void fill_left(const int* ptr, int n) {
        for (int k=0; k<n; ++k, ++ptr) {
            encountered[*ptr]=false;
        }
        return;
    }

    bool intersect_right(const int* ptr, int n, int i) {
        for (int k=0; k<n; ++k, ++ptr) {
            auto it=encountered.find(*ptr);
            if (it!=encountered.end() && !(it->second)) {
                it->second=true; // mark this subject as having been recorded already.
                qout.push_back(i+1); // get back to 1-based indexing.
                sout.push_back(*ptr);

                if (skip) { return true; } // quit if we only need a single hit for 'i'.
            }
        }
        return false;
    }

    Rcpp::List yield() const {
        return Rcpp::List::create(
            Rcpp::IntegerVector(qout.begin(), qout.end()),
            Rcpp::IntegerVector(sout.begin(), sout.end())
        );
    }

private:
    std::deque<int> qout, sout;
    std::unordered_map<int, bool> encountered;
    bool skip;
};

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
    Rcpp::IntegerVector subject_indices_right, Rcpp::IntegerVector subject_order_right,
    bool quit_on_first) 
{
    index_organizer left(query_hits_left, subject_hits_left, subject_indices_left, subject_order_left);
    index_organizer right(query_hits_right, subject_hits_right, subject_indices_right, subject_order_right);
    hits_collector output(quit_on_first);

    // Iterating across all queries.
    if (query_indices_left.size()!=query_indices_right.size()) {
        throw std::runtime_error("query index vectors should be the same length");
    }
    auto qlIt=query_indices_left.begin(), qrIt=query_indices_right.begin();
    const size_t N=query_indices_left.size();

    for (size_t i=0; i<N; ++i, ++qlIt, ++qrIt) {
        auto lqIt=left.search_query_region(*qlIt);
        if (!left.has_query_region(lqIt)) {
            continue;
        }
        auto rqIt=right.search_query_region(*qrIt);
        if (!right.has_query_region(rqIt)) {
            continue;
        }

        // We identify all subject regions that overlap with the left region of the current query.
        const auto& lopt=lqIt->second;
        const int* lptr=lopt.first;
        const int ln=lopt.second;

        for (int j=0; j<ln; ++j, ++lptr) {
            // Second search identifies all subject interactions with
            // a left region among the set of overlapped subject regions.
            auto lsIt=left.search_subject_region(*lptr);
            if (left.has_subject_region(lsIt)) {
                const auto& opt=lsIt->second;
                output.fill_left(opt.first, opt.second);
            }
        }

        if (output.no_left_hits()) {
            continue;
        }

        // We next identify all subject regions that overlap with the right region of the current query.
        const auto& ropt=rqIt->second;
        const int* rptr=ropt.first;
        const int rn=ropt.second;

        for (int j=0; j<rn; ++j, ++rptr) {
            // Second search identifies all subject interactions with
            // a right region among the overlapped subject regions. We then
            // intersect this with those we found for the left overlaps.
            auto rsIt=right.search_subject_region(*rptr);
            if (right.has_subject_region(rsIt)) { 
                const auto& opt=rsIt->second;
                if (output.intersect_right(opt.first, opt.second, i)) { break; } 
            }
        }

        output.reset_left_store();
    }

    return output.yield();
}
