#! /bin/bash

##
##  A function to extract a date stamp of type:
##+ '2013 December 19' out og a FORTRAN IERS program source file. In these
##+ files, there is a comment block (at header) of type:
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
##  Some of the source code, has the date revisions stamps in another format,
##+ i.e. '2015 16 June'. So, if we get an empty string from the first try, we
##+ will re-try getting the date string using this, alternative format.
##
##  Yet another format, is: 
##+ ' * Beth Stetzler, 31 May 2013, Updated reference and added 10 January[...]'
##+ so, in case of failure, we also check this one!
##
function extract_for_revision_date()
{
    ## first try, format is:
    ## '*  2013 December 19 G. Petit      Updated[...]'
    stamp=$(cat ${1} | \
        sed -n '/*  Revisions:/,/-------/p' | \
        grep '\*  [0-9]\{4\} \w*\s\{1,5\}[0-9]\{1,2\}' | \
        tail -1 | \
        grep -o '[0-9]\{4\} \w*\s\{1,5\}[0-9]\{1,2\}')
    ## second try, format is
    ## '*  2015 16 June     N. Stamatakos at[...]'
    if  test -z "$stamp" ; then
        stamp=$(cat ${1} | \
            sed -n '/*  Revisions:/,/-------/p' | \
            grep '\*  [0-9]\{4\} [0-9]\{1,2\} \w*' | \
            tail -1 | \
            grep -o '[0-9]\{4\} [0-9]\{1,2\} \w*')
    fi
    ## third try, format is
    ## ' * Beth Stetzler, 31 May 2013, Updated reference[...]'
    if  test -z "$stamp" ; then
        stamp=$(cat ${1} | \
            sed -n '/*  Revisions:/,/-------/p' | \
            grep '\* \w* \w*, [0-9]\{1,2\} \w* [0-9]\{4\},' | \
            tail -1 | \
            grep -o '[0-9]\{1,2\} \w* [0-9]\{4\}')
    fi
    if  test -z "$stamp" ; then
        echo "[ERROR] Failed to get date stamp from file \"${1}\"."
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
        grep '\* @version [0-9]\{2\}\.[0-9]\{2\}\.[0-9]\{4\}' | \
        tail -1 | \
        grep -o '[0-9]\{2\}\.[0-9]\{2\}\.[0-9]\{4\}')
    if  test -z "$stamp" ; then
        echo "[ERROR] Failed to get date stamp from file \"${1}\"."
        exit 1
    fi
    echo $stamp
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
echo "2. Extracting last revision date stamps"
for p in ${fprogs[@]} ; do
    f=$(basename $p)
    ## fortran
    if ! test -f $f ; then
        echo "  [${counter}/${#fprogs[@]}] File \"${f}\" does not exist!"
        exit 1
    else
        stamp=$(extract_for_revision_date ${f})
        echo "  [${counter}/${#fprogs[@]}] Last revision for ${f} is at \"$stamp\""
    fi
    ## cpp
    cppf=$(echo ${f,,} | sed 's/\.f/\.cpp/g')
    if ! test -f ${CPP_SRC_DIR}/${cppf} ; then
        echo "  [${counter}/${#fprogs[@]}] File \"${cppf}\" does not exist!"
        # exit 1
    else
        stamp=$(extract_cpp_revision_date ${CPP_SRC_DIR}/${cppf})
        echo "  [${counter}/${#fprogs[@]}] Last revision for ${cppf} is at \"$stamp\""
    fi
    let counter=counter+1
done

exit 0
