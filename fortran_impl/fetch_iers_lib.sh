#! /bin/bash

declare -A months=( ["jan"]="01" ["feb"]="02" ["mar"]="03" ["apr"]="04" \
["may"]="05" ["jun"]="06" ["jul"]="07" ["aug"]="08" ["sep"]="09" ["oct"]="10" \
["nov"]="11" ["dec"]="12")

FORTRAN_DIR=""
FORCE_DOWNLOAD=""

##
##  A function to extract a date stamp of type:
##+ '2013 December 19' out of a FORTRAN IERS program source file.
##
##  In these files, there is a comment block (at header) of type:
##  * [...]
##  *  Revisions:
##  *  2007 August   27 S. Lambert    Original code
##  *  2008 November 19 B.E.Stetzler  Added header and copyright
##  *  2008 November 20 B.E.Stetzler  provided a test case 
##  *                                 test case results
##  * [...]
##  *  2013 December 19 G. Petit      Updated parameter N to 30, 
##  *                                 updated table to 2013
##  *                                 test case results NOT UPDATED, TO BE DONE
##  *-----------------------------------------------------------------------
##+ so this function will extract the last revision date (in this case the
##+ string '2013 December 19'
##
##  Some of the source code, has the date revision stamps in another format,
##+ i.e. '2015 16 June'. So, if we get an empty string from the first try, we
##+ will re-try getting the date string using this, alternative format.
##
##  Yet another format, is: 
##+ ' * Beth Stetzler, 31 May 2013, Updated reference and added 10 January[...]'
##+ so, in case of failure, we also check this one!
##
##  In any case, the catpured date string will be returned in the format:
##+ '2013 December 19'
##
function extract_for_revision_date()
{
    if ! test -f "${1}" ; then
        echo "[ERROR] Failed to get date stamp from file \"${1}\"." 1>&2
        echo "[ERROR] File does not exist!" 1>&2
        exit 1
    fi

    ## first try, format is:
    ## '*  2013 December 19 G. Petit      Updated[...]'
    stamp=$(cat ${1} | \
        sed -n '/*  Revisions:/,/-------/p' | \
        grep '\*  [0-9]\{4\} \w*\s\{1,5\}[0-9]\{1,2\}' | \
        tail -1 | \
        grep -o '[0-9]\{4\} \w*\s\{1,5\}[0-9]\{1,2\}' 2>/dev/null)
    ## second try, format is
    ## '*  2015 16 June     N. Stamatakos at[...]'
    ## Output date stamp as '2015 June 16'
    if  test -z "$stamp" ; then
        stamp=$(cat ${1} | \
            sed -n '/*  Revisions:/,/-------/p' | \
            grep '\*  [0-9]\{4\} [0-9]\{1,2\} \w*' | \
            tail -1 | \
            sed -n 's/\*  \([0-9]\{4\}\) \([0-9]\{1,2\}\) \(\w*\).*/\1 \3 \2/p'\
            2>/dev/null)
    fi
    ## third try, format is
    ## ' * Beth Stetzler, 31 May 2013, Updated reference[...]'
    ## Output date stamp as '2013 May 31'
    if  test -z "$stamp" ; then
        stamp=$(cat ${1} | \
            sed -n '/*  Revisions:/,/-------/p' | \
            grep '\* \w* \w*, [0-9]\{1,2\} \w* [0-9]\{4\},' | \
            tail -1 | \
            sed -n 's/[^0-9]*\([0-9]\{1,2\}\) \(\w*\) \([0-9]\{4\}\).*/\3 \2 \1/p' \
            2>/dev/null)
    fi
    ## forth try is
    ## '*  This revision:  2014 November 7'
    ## This is used by CAL2JD and DAT (actually the SOFA sw)
    if  test -z "$stamp" ; then
        stamp=$(cat ${1} | \
            grep "This revision:" | \
            sed -n 's/^\*\s*This revision:\s*\([0-9]\{4\}\) \([a-z,A-Z]*\) \([0-9]\{1,2\}\)/\1 \2 \3/p'
            2>/dev/null)
    fi
    if  test -z "$stamp" ; then
        echo "[ERROR] Failed to get date stamp from file \"${1}\"." 1>&2
        exit 1
    fi
    echo $stamp
}

##
##  Extract the version date stamp from a c++ iers source code file.
##
function extract_cpp_revision_date()
{
    stamp=$(cat ${1} | \
        grep -o '\* @version\s*[0-9]\{2\}\.[0-9]\{2\}\.[0-9]\{4\}' | \
        tail -1 | \
        grep -o '[0-9]\{2\}\.[0-9]\{2\}\.[0-9]\{4\}')
    if  test -z "$stamp" ; then
        echo "[ERROR] Failed to get date stamp from file \"${1}\"." 1>&2
        exit 1
    fi
    echo $stamp
}

##
##  Compare two date strings of type:
##+ '2013 December 19' vs '19.12.2013'
##  If the dates match, 0 is returned, else 1.
##
compare_date_str()
{
    if ! test "$#" -eq 4 ; then
        # for i in "$*" ; do echo "\"$i\""; done
        echo "[ERROR] Invalid date strings to compare!"
        exit 1
    fi

    local mon=${2,,}
    if test "${#3}" -lt 2 ; then
        local day="0$3"
    else
        local day=$3
    fi
    local for_str=${day}.${months["${mon:0:3}"]}.${1}
    
    local cpp_str=${4}
    
    if test "$for_str" != "$cpp_str" ; then
        echo "[WARNING] Dates do not match (\"$for_str\" vs \"$cpp_str\")."
        return 1
    fi

    return 0
}

IERS_FTP_DIR="ftp://maia.usno.navy.mil/conventions/2010/2010_update"
CPP_SRC_DIR="../src"

fprogs=("chapter5/software/FCNNUT.F" \
"chapter5/software/FUNDARG.F" \
"chapter5/software/PMSDNUT2.F" \
"chapter5/software/UTLIBR.F" \
"chapter7/software/IERS_CMP_2015.F" \
"chapter7/software/ARG2.F" \
"chapter8/software/CNMTX.F" \
"chapter8/software/ORTHO_EOP.F" \
"chapter8/software/RG_ZONT2.F" \
"chapter9/software/APG.F" \
"chapter9/software/FCUL_A.F" \
"chapter9/software/FCUL_B.F" \
"chapter9/software/FCUL_ZD_HPA.F" \
"chapter9/software/GMF.F" \
"chapter9/software/GPT.F" \
"chapter9/software/GPT2.F" \
"chapter9/software/VMF1.F" \
"chapter9/software/VMF1_HT.F")

detide_progs=("chapter7/software/dehanttideinel/CAL2JD.F" \
"chapter7/software/dehanttideinel/DAT.F" \
"chapter7/software/dehanttideinel/DEHANTTIDEINEL.F" \
"chapter7/software/dehanttideinel/NORM8.F" \
"chapter7/software/dehanttideinel/SPROD.F" \
"chapter7/software/dehanttideinel/ST1IDIU.F" \
"chapter7/software/dehanttideinel/ST1ISEM.F" \
"chapter7/software/dehanttideinel/ST1L1.F" \
"chapter7/software/dehanttideinel/STEP2DIU.F" \
"chapter7/software/dehanttideinel/STEP2LON.F" \
"chapter7/software/dehanttideinel/ZERO_VEC8.F")

fprogs+=( "${detide_progs[@]}" )
echo "1. Downloading IERS source code from ${IERS_FTP_DIR}"
counter=1
for p in ${fprogs[@]} ; do
    f=$(basename $p)
    if ! test -f ${f} ; then
        echo "  [${counter}/${#fprogs[@]}] Downloading ${IERS_FTP_DIR}/${p} ...."
        wget -q -O ${f} ${IERS_FTP_DIR}/${p} || {\
            rm -f ${f}
            echo "[ERROR] Failed to download file \"${IERS_FTP_DIR}/${p}\""
            exit 1
        }
    else
        echo "  [${counter}/${#fprogs[@]}] File \"${f}\" already exists; skipping download"
    fi
    let counter=counter+1
done

counter=1
echo "2. Extracting and comparing last revision date stamps"
set -e
for p in ${fprogs[@]} ; do
    f=$(basename $p)
    # echo "Comparing $f"
    ## fortran
    if ! test -f $f ; then
        echo "  [${counter}/${#fprogs[@]}] File \"${f}\" does not exist!"
        exit 1
    else
        for_stamp=$(extract_for_revision_date ${f})
        # echo "  [${counter}/${#fprogs[@]}] Last revision for ${f} is at \"$for_stamp\""
    fi
    ## cpp
    cppf=$(echo ${f,,} | sed 's/\.f/\.cpp/g')
    if ! test -f ${CPP_SRC_DIR}/${cppf} ; then
        echo "  [${counter}/${#fprogs[@]}] File \"${cppf}\" does not exist!"
        # exit 1
    else
        cpp_stamp=$(extract_cpp_revision_date ${CPP_SRC_DIR}/${cppf})
        # echo "  [${counter}/${#fprogs[@]}] Last revision for ${cppf} is at \"$cpp_stamp\""
        compare_date_str $for_stamp $cpp_stamp || echo "[WARNING] Program name: ${cppf/.cppp/}"
    fi
    let counter=counter+1
done
set +e

echo "3. Editing FORTRAN programs and needed data"
##  In DAT.F the call to CAL2JD is: 'call iau_CAL2JD (...)'; this will not work
##+ it needs to be changed to: 'call CAL2JD (...)'
sed -i 's/call iau_CAL2JD/call CAL2JD/g' DAT.F
##  Link the 'data/gpt2_5.grd' in the FORTRAN folder
ln -s 

echo "4. Compiling FORTRAN programs"
make || (echo "[ERROR] Failed to compile FORTRAN source code." && exit 1)

exit 0
