################################################################################
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
#
# Author:  Douglas Slotta
#
################################################################################

flex(vcf_scanner vcf-scan.ll ${CMAKE_CURRENT_BINARY_DIR}/vcf-scan.cpp FLAGS "-Cf")
#SET_SOURCE_FILES_PROPERTIES(${FLEX_vcf_scanner_OUTPUTS})
#add_flex_bison_dependency(vcf_scanner)

set (vcf2scidb_src
  vcf2scidb.cpp
  scidb-writers.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/vcf-scan.cpp
  )

file(GLOB vcf2scidb_inc "*.hpp" "*.ll")
#set(vcf2scidb_inc scidb-writers.hpp)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(vcf2scidb ${vcf2scidb_src} ${vcf2scidb_inc})
extractDebugInfo("${GENERAL_OUTPUT_DIRECTORY}" "vcf2scidb" vcf2scidb)
set_target_properties(vcf2scidb PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${GENERAL_OUTPUT_DIRECTORY})
target_link_libraries(vcf2scidb
    ${Boost_LIBRARIES}
)

set_target_properties(vcf2scidb
    PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE
    CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE
    INSTALL_RPATH ${DEFAULT_RPATH}
    COMPILE_FLAGS "-std=c++0x"
)
