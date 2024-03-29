# Require dot, treat the other components as optional
find_package(Doxygen REQUIRED dot)

# equals doxyfile_common.in
set(DOXYGEN_HTML_HEADER            "${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/header.html")
set(DOXYGEN_HTML_FOOTER            "${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/footer.html")
set(DOXYGEN_GENERATE_TREEVIEW      "YES")
set(DOXYGEN_DISABLE_INDEX          "NO")
set(DOXYGEN_FULL_SIDEBAR           "NO")
set(DOXYGEN_HTML_EXTRA_STYLESHEET  "${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/doxygen-awesome.css")
set(DOXYGEN_HTML_EXTRA_FILES       "${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/doxygen-awesome-darkmode-toggle.js;"
                                   "${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/doxygen-awesome-interactive-toc.js")
set(DOXYGEN_GENERATE_LATEX         "NO")
set(DOXYGEN_EXTRACT_PRIVATE        "YES")
set(DOXYGEN_USE_MATHJAX            "YES")
set(DOXYGEN_MATHJAX_RELPATH        "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/")
set(DOXYGEN_HAVE_DOT               "YES")
set(DOXYGEN_COLLABORATION_GRAPH    "NO")
set(DOXYGEN_DOT_IMAGE_FORMAT       "svg")
set(DOXYGEN_MAX_DOT_GRAPH_DEPTH    "2")
set(DOXYGEN_DOT_MULTI_TARGETS      "YES")
set(DOXYGEN_DOT_TRANSPARENT        "YES")
set(DOXYGEN_WARN_IF_UNDOCUMENTED   "NO")
set(DOXYGEN_WARN_IF_INCOMPLETE_DOC "NO")
set(DOXYGEN_WARN_NO_PARAMDOC       "NO")

# equals doxygen.in
set(DOXYGEN_ PROJECT_BRIEF         "Fast Matrix/Vector Library")
# INPUT          set automatically
# PROJECT_NAME   set automatically
# PROJECT_NUMBER set automatically
set(DOXYGEN_ GENERATE_TAGFILE      "${CMAKE_CURRENT_BINARY_DIR}/html/fmatvec.tag")
# TAGFILES not needed

doxygen_add_docs(doc
  ${PROJECT_SOURCE_DIR}/fmatvec
  COMMENT "Generate doxygen doc"
)

# Check for node (implying mathjax-node is installed) and add check as POST_BUILD
# command to target doc
find_program(NODE node)
if(NODE)
  add_custom_command(TARGET doc
                     POST_BUILD
                     COMMAND ${NODE} ${CMAKE_CURRENT_SOURCE_DIR}/checkHtml/checkMathJax.js ${CMAKE_CURRENT_BINARY_DIR}
  )
else()
  message(WARNING "mathjax-node/node not found: skipping check of math formula generation for doxygen output")
endif()

# Configure doxyfile_common
set(prefix ${CMAKE_INSTALL_PREFIX})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/doxyfile_common.in ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/doxyfile_common)

# copy stuff that is needed locally
add_custom_target(doc_copy
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/header.html ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/footer.html ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/doxygen-awesome.css ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/doxygen-awesome-darkmode-toggle.js ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/doxygen-awesome-interactive-toc.js ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/doxy_template/doxy.mk     ${CMAKE_CURRENT_BINARY_DIR}/doxy_template/
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
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/checkHtml/
        DESTINATION share/fmatvec/checkHtml
        EXCLUDE_FROM_ALL
        COMPONENT doc
        PATTERN "*.sh" PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
        PATTERN "Makefile.*" EXCLUDE
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
