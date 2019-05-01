#include "Rcpp.h"

#include <deque>
#include <algorithm>
#include <stdexcept>
#include <string>

void flip(Rcpp::String& left_chr, Rcpp::String& right_chr,
    int& left_start, int& right_start,
    int& left_end, int& right_end) 
{
    if (left_chr > right_chr) {
        std::swap(left_chr, right_chr);
        std::swap(left_start, right_start);
        std::swap(left_end, right_end);

    } else if (left_chr==right_chr) {
        if (left_start > right_start) {
            std::swap(left_start, right_start);
            std::swap(left_end, right_end);

        } else if (left_start==right_start) {
            if (left_end > right_end) {
                std::swap(left_end, right_end);
            }
        }
    }
    return;
}

// [[Rcpp::export(rng=false)]]
Rcpp::List bounding_box(Rcpp::IntegerVector runs, Rcpp::StringVector values,
    Rcpp::IntegerVector index1, Rcpp::StringVector ref_chr1, Rcpp::IntegerVector ref_start1, Rcpp::IntegerVector ref_end1,
    Rcpp::IntegerVector index2, Rcpp::StringVector ref_chr2, Rcpp::IntegerVector ref_start2, Rcpp::IntegerVector ref_end2,
    Rcpp::LogicalVector reflect)
{
    BEGIN_RCPP
    const size_t ngroups=runs.size();
    if (runs.size()!=values.size()) {
        throw std::runtime_error("'runs' and 'values' should have the same length");
    }
    if (index1.size()!=index2.size()) {
        throw std::runtime_error("'index1' and 'index2' should have the same length");
    }
    if (reflect.size()!=1) {
        throw std::runtime_error("'reflect' should be a logical scalar");
    }
    const bool do_reflect=reflect[0];

    // Setting up output vectors.
    Rcpp::StringVector cout_left(ngroups), cout_right(ngroups);
    Rcpp::IntegerVector sout_left(ngroups), eout_left(ngroups), sout_right(ngroups), eout_right(ngroups);

    auto i1it=index1.begin();
    auto i2it=index2.begin();
    for (size_t r=0; r<ngroups; ++r) {
        const int currun=runs[r];
        if (currun==0) { throw std::runtime_error("empty group"); }

        // Setting up starting values from first element in the group.
        const int a1=*i1it-1, a2=*i2it-1;
        Rcpp::String left_chr=ref_chr1[a1], right_chr=ref_chr2[a2];
        int left_start=ref_start1[a1], left_end=ref_end1[a1],
            right_start=ref_start2[a2], right_end=ref_end2[a2];

        if (do_reflect) {
            flip(left_chr, right_chr, left_start, right_start, left_end, right_end);
        }
        ++i1it;
        ++i2it;

        // Running through the remaining elements and doing a min/max on the starts and ends.
        for (int i=1; i<currun; ++i, ++i1it, ++i2it) {
            const int a1=*i1it-1, a2=*i2it-1;
            Rcpp::String chr1=ref_chr1[a1], chr2=ref_chr2[a2];
            int left_start_candidate=ref_start1[a1];
            int left_end_candidate=ref_end1[a1];
            int right_start_candidate=ref_start2[a2];
            int right_end_candidate=ref_end2[a2];

            if (do_reflect) {
                flip(chr1, chr2, left_start_candidate, right_start_candidate,
                    left_end_candidate, right_end_candidate);
            }

            if (chr1!=left_chr || chr2!=right_chr) {
                throw std::runtime_error(std::string("multiple chromosomes for group '")
                    + Rcpp::as<std::string>(values[r]) + "'");
            }

            left_start=std::min(left_start, left_start_candidate);
            right_start=std::min(right_start, right_start_candidate);
            left_end=std::max(left_end, left_end_candidate);
            right_end=std::max(right_end, right_end_candidate);
        }

        cout_left[r]=left_chr;
        cout_right[r]=right_chr;
        sout_left[r]=left_start;
        eout_left[r]=left_end;
        sout_right[r]=right_start;
        eout_right[r]=right_end;
    }
    
    return Rcpp::List::create(
        Rcpp::List::create(cout_left, sout_left, eout_left),
        Rcpp::List::create(cout_right, sout_right, eout_right)
    );
    END_RCPP
}
