#!/bin/bash
# Stage an input file for metaMS (rule metab_metams_one): usage metams_stage_mzml.sh <in> <out>
#
# Thermo files (.raw, or mzML from ThermoRawFileParser) are copied as is. Other mzML (msconvert
# output from Agilent, Bruker, Sciex, ... or mzXML) is rewritten, because CoreMS's mzML reader
# assumes Thermo conventions:
#   - spectrum ids become Thermo-style "controllerType=0 controllerNumber=1 scan=N", N = index + 1
#     (CoreMS drops those id fields unconditionally and fails when they are missing); precursor
#     spectrumRefs are remapped to match
#   - each spectrum gets a "filter string" (MS:1000512) of "ms1"/"ms2" from its ms level; metaMS
#     selects MS2 scans by matching the scan translator's scan_filter against it (Thermo's filter
#     string otherwise; config/metams/generic_no_filter.toml matches "ms2")
#   - the indexedmzML wrapper is dropped, since its byte offsets no longer match; the mzML
#     element gets the default namespace that pymzml needs if it lacks one
set -euo pipefail
in=$1; out=$2
if [[ ! "${in,,}" == *.mzml ]] || grep -q -m1 'controllerType=' "$in"; then
    cp "$in" "$out"
    exit 0
fi
{
    echo '<?xml version="1.0" encoding="utf-8"?>'
    sed -n '/<mzML /,/<\/mzML>/p' "$in" | awk '
        # default namespace, unless the file already declares it
        !ns && /<mzML / { if ($0 !~ /<mzML [^>]*xmlns="/) sub(/<mzML /, "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" "); ns = 1 }
        # spectrum ids -> sequential Thermo-style scan numbers (index + 1), whatever the
        # vendor id was (msconvert "scan=N", Agilent "scanId=N", ...); precursor
        # spectrumRefs point at earlier spectra, so one pass can remap them
        /<spectrum / && match($0, / id="[^"]*"/) {
            old = substr($0, RSTART + 5, RLENGTH - 6)
            match($0, / index="[0-9]+"/); n = substr($0, RSTART + 8, RLENGTH - 9) + 1
            new = "controllerType=0 controllerNumber=1 scan=" n; ids[old] = new
            sub(/ id="[^"]*"/, " id=\"" new "\"")
        }
        /spectrumRef="/ && match($0, /spectrumRef="[^"]*"/) {
            ref = substr($0, RSTART + 13, RLENGTH - 14)
            if (ref in ids) sub(/spectrumRef="[^"]*"/, "spectrumRef=\"" ids[ref] "\"")
        }
        { print }
        # filter string "ms<level>" after each ms level cvParam
        /accession="MS:1000511"/ {
            match($0, /^ */); ind = substr($0, 1, RLENGTH)
            match($0, /value="[0-9]+"/); lvl = substr($0, RSTART + 7, RLENGTH - 8)
            print ind "<cvParam cvRef=\"MS\" accession=\"MS:1000512\" name=\"filter string\" value=\"ms" lvl "\"/>"
        }'
} > "$out"
