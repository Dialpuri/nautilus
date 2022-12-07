/*! \file nautilus-target.h nautilus library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#ifndef NAUTILUS_TARGET_H
#define NAUTILUS_TARGET_H


#include "nucleicacid_db.h"
#include "nautilus-ss-find.h"


class NucleicAcidTarget {
 public:
  typedef std::pair<clipper::Coord_orth,clipper::Coord_orth> Coords;
  typedef std::vector<Coords> Target;
  typedef std::vector<clipper::Coord_orth> Standard;
  NucleicAcidTarget() {}
  void init( const float c_hi[][3], const float c_lo[][3], const float c_repr[3][3], const int ncoord );
  void init_stats( const clipper::Xmap<float>& xmap );
  const Target& target() const { return target_; }
  const Standard& standard() const { return standard_; }
  double radius() const;

  float score_min( const clipper::Xmap<float>& xmap, const clipper::RTop_orth& rtop ) const;
  float score_sum( const clipper::Xmap<float>& xmap, const clipper::RTop_orth& rtop ) const;
  float cutoff_min( double p ) const;
  float cutoff_sum( double p ) const;

 private:
  Target target_;
  std::vector<clipper::Coord_orth> standard_;
  std::vector<float> smin, ssum;
};

class NucleicAcidTargets {
 public:
  NucleicAcidTargets();
  void add_pdb( const clipper::String& file );
  void init_stats( const clipper::Xmap<float>& xmap );
  const NucleicAcidDB::Chain& db() const { return nadb; } 
  static void superpose_sugar( NucleicAcidDB::Chain& frag, int posn, const NucleicAcidDB::NucleicAcid& na );
  float score_sugar( const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na ) const;
  float score_phosphate( const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na1, const NucleicAcidDB::NucleicAcid& na2 ) const;
  NucleicAcidDB::NucleicAcid next_na_group( const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na ) const;
  NucleicAcidDB::NucleicAcid prev_na_group( const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na ) const;
  const NucleicAcidTarget& target_sugar() const { return target_s; }
  const NucleicAcidTarget& target_phosphate() const { return target_p; }

  const NucleicAcidDB::Chain join_sugars( const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na1, const NucleicAcidDB::NucleicAcid& na2, int len, double rmsdlim ) const;

  const clipper::MiniMol phosphate( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol, const clipper::MiniMol& mol_pho );
  const clipper::MiniMol find( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol, int nsugar, int nphosp, double step );
  const clipper::MiniMol grow( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol, int ngrow, double fcut ) const;
  const clipper::MiniMol link( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol ) const;
  const clipper::MiniMol prune( const clipper::MiniMol& mol ) const;
  const clipper::MiniMol rebuild_chain( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol ) const;

//  Use correlation based function to score the fragment
  const float calculate_correlation(const clipper::Xmap<float>& xmap, const clipper::Xmap<float>& mask) const;
  const float score_na_fragment(const clipper::Xmap<float>& xmap, NucleicAcidDB::NucleicAcid& fragment) const;
  const float score_binucleotide(const clipper::Xmap<float>& xmap, NucleicAcidDB::Chain& current_fragment );


  const void dump_ed_around_fragment(const clipper::Xmap<float>& xmap, NucleicAcidDB::NucleicAcid fragment);
  void dump_search_atoms(std::vector<std::pair<std::string, clipper::Coord_orth>> input_atoms, std::basic_string<char> path, std::basic_string<char> name);

    template <typename T>
    std::string to_string_with_precision(const T a_value, const int n = 4)
    {
        std::ostringstream out;
        out.precision(n);
        out << std::fixed << a_value;
        std::string out_str(6, '\0');

        for (int i = 0; i < out.str().size(); i++) {
            out_str[i] += out.str()[i];
        }
        return out_str;
    }

 private:
  NucleicAcidDB::Chain nadb;
  NucleicAcidDB::NucleicAcid narepr;
  NucleicAcidTarget target_s, target_p;
  std::vector<SearchResult> found_s, found_p;


};

#endif
