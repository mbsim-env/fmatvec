# Require dot, treat the other components as optional
find_package(Doxygen
             REQUIRED dot)

set(DOXYGEN_PROJECT_BRIEF "Fast Matrix/Vector Library")
set(DOXYGEN_GENERATE_LATEX      "NO")
set(DOXYGEN_EXTRACT_PRIVATE     "YES")
set(DOXYGEN_USE_MATHJAX         "YES")
set(DOXYGEN_HAVE_DOT            "YES")
set(DOXYGEN_COLLABORATION_GRAPH "NO")
set(DOXYGEN_DOT_IMAGE_FORMAT    "svg")
set(DOXYGEN_MAX_DOT_GRAPH_DEPTH "2")
set(DOXYGEN_DOT_MULTI_TARGETS   "YES")
set(DOXYGEN_DOT_TRANSPARENT     "YES")
set(DOXYGEN_GENERATE_TAGFILE    "${CMAKE_CURRENT_BINARY_DIR}/html/fmatvec.tag")

set(DOXYGEN_HTML_HEADER           "${CMAKE_CURRENT_BINARY_DIR}/doxy_template/header.html")
set(DOXYGEN_HTML_FOOTER           "${CMAKE_CURRENT_BINARY_DIR}/doxy_template/footer.html")
set(DOXYGEN_HTML_EXTRA_STYLESHEET "${CMAKE_CURRENT_BINARY_DIR}/doxy_template/customdoxygen.css")

doxygen_add_docs(doc
    ${PROJECT_SOURCE_DIR}/fmatvec
    COMMENT "Generate doxygen doc"
)

# copy stuff that is needed locally
add_custom_target(doc_copy
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/customdoxygen.css ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/footer.html       ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/header.html       ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/doxy.mk           ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
    )

add_dependencies(doc doc_copy)

# make doc/clean automatically provided via package DOXYGEN

# unfortunately, this install does not work but even breaks the build => TODO: rework
# # make install for html doc
# install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/          DESTINATION share/doc/fmatvec/)
# install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/ DESTINATION share/fmatvec/doxy_template)
# 
