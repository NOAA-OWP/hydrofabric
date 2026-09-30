#' hydrofabric: Build Hydrologic Modeling Fabrics from Reference Networks
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats na.omit setNames
#' @importFrom utils modifyList read.csv write.csv
#' @importFrom data.table :=
#' @importFrom dplyr all_of
NULL

# Non-standard-evaluation column names used across dplyr/data.table pipelines.
# Declared here to silence R CMD check 'no visible binding' NOTEs without
# qualifying every column reference in the codebase.
utils::globalVariables(c(
  ".", ".area", ".data", ".feature_id", ".I", ".N", ".part_id", ".SD",
  "absorbed_recs", "aggregated_id", "areasqkm", "attr_source", "base", "becomes", "bf_area", "BtmWdth", "Bw",
  "child_area", "ChSlp", "comid", "COMID", "divide_id", "downstream_id", "downstream_x", "downstream_y",
  "ds_id", "e_end", "e_start", "end", "event_identifier", "event_REACH_meas", "event_REACHCODE", "event_split_fID",
  "fac_value", "feature_id", "fID", "final_divide_id", "final_flowpath_id", "fl_id", "fl1", "fl2",
  "flowline_hydroseq", "flowline_id", "flowline_incr", "flowline_lengthkm", "flowline_toid", "flowpath_id", "flowpath_toid", "fp_idx",
  "from_fp_int", "geom_len", "geometry", "group_id", "has_divide", "hf_id", "hf_id_part", "hl_class",
  "hl_dn", "hl_reference", "hl_un", "hl2", "hydroseq", "i.ID", "i.new_lp", "ID",
  "ID_Hydroseq", "ID_new", "incremental_area_id", "ind", "is_rep", "is_terminal", "join_category", "joined_from_flowpath_id",
  "joined_to_flowpath_id", "keep", "lake_id", "length_m", "lengthkm", "levelpathid", "mainstem", "mC1",
  "mC2", "member", "member_comid", "member_flowpath_id", "min_hs", "MusK", "MusX", "N",
  "nCC", "new_becomes", "new_fp_int", "new_int", "nhd_comid", "nP", "old_int", "orig_flowpath_id",
  "orig_levelpathid", "owp_roughness_bathy", "owp_roughness_overbank",
  "owp_tw_bf", "owp_y_bf", "parent_div", "parent_divide_id", "part", "part_num",
  "piece", "pieces", "poi_id", "REACH_meas", "rec_id", "rec_id_list", "rec_toid", "recipient",
  "reconciled_id", "reconciled_toid", "refactored_flowpath_id", "refactored_id", "reference_id", "seg_id", "sent", "set", "slope",
  "So", "split_part", "start", "to_becomes", "to_fp_int", "toCOMID", "toid", "toID",
  "toID_new", "TopWdth", "TopWdthCC", "toset", "tot_drainage_area", "tot_drainage_areasqkm", "total_dasqkm", "toVPUID",
  "Tw", "type", "upstream_count", "vpuid", "VPUID", "wbareacomi", "Y"
))
