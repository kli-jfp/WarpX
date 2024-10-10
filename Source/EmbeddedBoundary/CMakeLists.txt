find_package(Boost 1.66.0 COMPONENTS headers)

foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(lib_${SD}
      PRIVATE
        Enabled.cpp
        WarpXInitEB.cpp
        WarpXFaceExtensions.cpp
        WarpXFaceInfoBox.H
        Enabled.cpp
        PolygonXYIF.H
    )
endforeach()
