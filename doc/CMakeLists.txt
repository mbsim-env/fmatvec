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

set(prefix ${CMAKE_INSTALL_PREFIX})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/doxyfile_common.in ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/doxyfile_common)

# copy stuff that is needed locally
add_custom_target(doc_copy
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/customdoxygen.css ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/footer.html       ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/header.html       ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/doxy.mk           ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/doxy-boot.js      ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
    COMMENT "Populating ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template"
    )

add_dependencies(doc doc_copy)

# The installation of the documentation must not take place in conjunction with main
# install. Therefore exclude it from ALL installation and make it a component "doc" install
# triggered via explicit cmake command "cmake -DCMAKE_INSTALL_COMPONENT=doc ..." as
# custom target.
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
        DESTINATION share/doc/fmatvec/
        EXCLUDE_FROM_ALL
        COMPONENT doc
        )
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
        DESTINATION share/fmatvec/doxy_template
        EXCLUDE_FROM_ALL
        COMPONENT doc
        )

# custom target doc-all (alias for "doc"
add_custom_target(doc-all DEPENDS doc)

# custom target doc-install
add_custom_target(doc-install
        DEPENDS doc-all
        COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_COMPONENT=doc -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
        COMMENT "Installing documentation to ${CMAKE_INSTALL_PREFIX}/share"
        )

# custom target doc-clean
add_custom_target(doc-clean
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/html
        COMMENT "Removing recursively ${CMAKE_CURRENT_BINARY_DIR}/html"
        )
