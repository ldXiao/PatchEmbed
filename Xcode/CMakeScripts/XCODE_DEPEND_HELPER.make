# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.algorithms.Debug:
PostBuild.bcclean.Debug:
/Users/vector_cat/gits/bcClean/Xcode/Debug/libbcclean.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/Debug/libbcclean.a


PostBuild.embree.Debug:
PostBuild.algorithms.Debug: /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libembree3.a
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libembree3.a:\
	/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/bcclean.build/Debug/algorithms.build/Objects-normal/libalgorithms.a
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libembree3.a


PostBuild.glad.Debug:
/Users/vector_cat/gits/bcClean/Xcode/glad/Debug/libglad.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/glad/Debug/libglad.a


PostBuild.glfw.Debug:
/Users/vector_cat/gits/bcClean/Xcode/glfw/src/Debug/libglfw3.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/glfw/src/Debug/libglfw3.a


PostBuild.imgui.Debug:
/Users/vector_cat/gits/bcClean/Xcode/imgui/Debug/libimgui.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/imgui/Debug/libimgui.a


PostBuild.lexers.Debug:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/liblexers.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/liblexers.a


PostBuild.math.Debug:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libmath.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libmath.a


PostBuild.predicates.Debug:
/Users/vector_cat/gits/bcClean/Xcode/predicates/Debug/libpredicates.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/predicates/Debug/libpredicates.a


PostBuild.simd.Debug:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libsimd.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libsimd.a


PostBuild.sys.Debug:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libsys.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libsys.a


PostBuild.tasking.Debug:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libtasking.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libtasking.a


PostBuild.test_bin.Debug:
PostBuild.bcclean.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.embed.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.json.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.igl_cgal.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.CGAL_Core.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.CGAL.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.igl_embree.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.embree.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.simd.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.lexers.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.sys.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.math.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.tasking.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.igl_opengl_glfw_imgui.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.igl_opengl_glfw.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.igl_opengl.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.imgui.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.glfw.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.glad.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.igl_predicates.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.igl.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.igl_common.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.predicates.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.cxxopts.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.nanoflann.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
PostBuild.pybind11.Debug: /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin
/Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin:\
	/Users/vector_cat/gits/bcClean/Xcode/Debug/libbcclean.a\
	/usr/local/lib/libboost_thread-mt.a\
	/usr/local/lib/libboost_system-mt.a\
	/usr/local/lib/libboost_chrono-mt.a\
	/usr/local/lib/libboost_date_time-mt.a\
	/usr/local/lib/libboost_atomic-mt.a\
	/usr/local/lib/libgmp.dylib\
	/usr/local/lib/libmpfr.dylib\
	/usr/local/lib/libboost_thread-mt.a\
	/usr/local/lib/libboost_chrono-mt.a\
	/usr/local/lib/libboost_date_time-mt.a\
	/usr/local/lib/libboost_atomic-mt.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libembree3.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libsimd.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/liblexers.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libsys.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libmath.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libtasking.a\
	/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/System/Library/Frameworks/OpenGL.framework/OpenGL.tbd\
	/Users/vector_cat/gits/bcClean/Xcode/imgui/Debug/libimgui.a\
	/Users/vector_cat/gits/bcClean/Xcode/glfw/src/Debug/libglfw3.a\
	/Users/vector_cat/gits/bcClean/Xcode/glad/Debug/libglad.a\
	/Users/vector_cat/gits/bcClean/Xcode/predicates/Debug/libpredicates.a\
	/usr/local/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7m.dylib
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/Debug/test_bin


PostBuild.algorithms.Release:
PostBuild.bcclean.Release:
/Users/vector_cat/gits/bcClean/Xcode/Release/libbcclean.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/Release/libbcclean.a


PostBuild.embree.Release:
PostBuild.algorithms.Release: /Users/vector_cat/gits/bcClean/Xcode/embree/Release/libembree3.a
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libembree3.a:\
	/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/bcclean.build/Release/algorithms.build/Objects-normal/libalgorithms.a
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Release/libembree3.a


PostBuild.glad.Release:
/Users/vector_cat/gits/bcClean/Xcode/glad/Release/libglad.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/glad/Release/libglad.a


PostBuild.glfw.Release:
/Users/vector_cat/gits/bcClean/Xcode/glfw/src/Release/libglfw3.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/glfw/src/Release/libglfw3.a


PostBuild.imgui.Release:
/Users/vector_cat/gits/bcClean/Xcode/imgui/Release/libimgui.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/imgui/Release/libimgui.a


PostBuild.lexers.Release:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/liblexers.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Release/liblexers.a


PostBuild.math.Release:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libmath.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Release/libmath.a


PostBuild.predicates.Release:
/Users/vector_cat/gits/bcClean/Xcode/predicates/Release/libpredicates.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/predicates/Release/libpredicates.a


PostBuild.simd.Release:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libsimd.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Release/libsimd.a


PostBuild.sys.Release:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libsys.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Release/libsys.a


PostBuild.tasking.Release:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libtasking.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Release/libtasking.a


PostBuild.test_bin.Release:
PostBuild.bcclean.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.embed.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.json.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.igl_cgal.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.CGAL_Core.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.CGAL.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.igl_embree.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.embree.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.simd.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.lexers.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.sys.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.math.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.tasking.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.igl_opengl_glfw_imgui.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.igl_opengl_glfw.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.igl_opengl.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.imgui.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.glfw.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.glad.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.igl_predicates.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.igl.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.igl_common.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.predicates.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.cxxopts.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.nanoflann.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
PostBuild.pybind11.Release: /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin
/Users/vector_cat/gits/bcClean/Xcode/Release/test_bin:\
	/Users/vector_cat/gits/bcClean/Xcode/Release/libbcclean.a\
	/usr/local/lib/libboost_thread-mt.a\
	/usr/local/lib/libboost_system-mt.a\
	/usr/local/lib/libboost_chrono-mt.a\
	/usr/local/lib/libboost_date_time-mt.a\
	/usr/local/lib/libboost_atomic-mt.a\
	/usr/local/lib/libgmp.dylib\
	/usr/local/lib/libmpfr.dylib\
	/usr/local/lib/libboost_thread-mt.a\
	/usr/local/lib/libboost_chrono-mt.a\
	/usr/local/lib/libboost_date_time-mt.a\
	/usr/local/lib/libboost_atomic-mt.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libembree3.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libsimd.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/Release/liblexers.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libsys.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libmath.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libtasking.a\
	/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/System/Library/Frameworks/OpenGL.framework/OpenGL.tbd\
	/Users/vector_cat/gits/bcClean/Xcode/imgui/Release/libimgui.a\
	/Users/vector_cat/gits/bcClean/Xcode/glfw/src/Release/libglfw3.a\
	/Users/vector_cat/gits/bcClean/Xcode/glad/Release/libglad.a\
	/Users/vector_cat/gits/bcClean/Xcode/predicates/Release/libpredicates.a\
	/usr/local/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7m.dylib
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/Release/test_bin


PostBuild.algorithms.MinSizeRel:
PostBuild.bcclean.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/libbcclean.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/libbcclean.a


PostBuild.embree.MinSizeRel:
PostBuild.algorithms.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libembree3.a
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libembree3.a:\
	/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/bcclean.build/MinSizeRel/algorithms.build/Objects-normal/libalgorithms.a
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libembree3.a


PostBuild.glad.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/glad/MinSizeRel/libglad.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/glad/MinSizeRel/libglad.a


PostBuild.glfw.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/glfw/src/MinSizeRel/libglfw3.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/glfw/src/MinSizeRel/libglfw3.a


PostBuild.imgui.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/imgui/MinSizeRel/libimgui.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/imgui/MinSizeRel/libimgui.a


PostBuild.lexers.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/liblexers.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/liblexers.a


PostBuild.math.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libmath.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libmath.a


PostBuild.predicates.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/predicates/MinSizeRel/libpredicates.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/predicates/MinSizeRel/libpredicates.a


PostBuild.simd.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libsimd.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libsimd.a


PostBuild.sys.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libsys.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libsys.a


PostBuild.tasking.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libtasking.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libtasking.a


PostBuild.test_bin.MinSizeRel:
PostBuild.bcclean.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.embed.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.json.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.igl_cgal.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.CGAL_Core.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.CGAL.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.igl_embree.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.embree.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.simd.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.lexers.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.sys.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.math.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.tasking.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.igl_opengl_glfw_imgui.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.igl_opengl_glfw.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.igl_opengl.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.imgui.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.glfw.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.glad.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.igl_predicates.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.igl.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.igl_common.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.predicates.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.cxxopts.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.nanoflann.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
PostBuild.pybind11.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin
/Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin:\
	/Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/libbcclean.a\
	/usr/local/lib/libboost_thread-mt.a\
	/usr/local/lib/libboost_system-mt.a\
	/usr/local/lib/libboost_chrono-mt.a\
	/usr/local/lib/libboost_date_time-mt.a\
	/usr/local/lib/libboost_atomic-mt.a\
	/usr/local/lib/libgmp.dylib\
	/usr/local/lib/libmpfr.dylib\
	/usr/local/lib/libboost_thread-mt.a\
	/usr/local/lib/libboost_chrono-mt.a\
	/usr/local/lib/libboost_date_time-mt.a\
	/usr/local/lib/libboost_atomic-mt.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libembree3.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libsimd.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/liblexers.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libsys.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libmath.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libtasking.a\
	/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/System/Library/Frameworks/OpenGL.framework/OpenGL.tbd\
	/Users/vector_cat/gits/bcClean/Xcode/imgui/MinSizeRel/libimgui.a\
	/Users/vector_cat/gits/bcClean/Xcode/glfw/src/MinSizeRel/libglfw3.a\
	/Users/vector_cat/gits/bcClean/Xcode/glad/MinSizeRel/libglad.a\
	/Users/vector_cat/gits/bcClean/Xcode/predicates/MinSizeRel/libpredicates.a\
	/usr/local/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7m.dylib
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/test_bin


PostBuild.algorithms.RelWithDebInfo:
PostBuild.bcclean.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/libbcclean.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/libbcclean.a


PostBuild.embree.RelWithDebInfo:
PostBuild.algorithms.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libembree3.a
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libembree3.a:\
	/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/bcclean.build/RelWithDebInfo/algorithms.build/Objects-normal/libalgorithms.a
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libembree3.a


PostBuild.glad.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/glad/RelWithDebInfo/libglad.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/glad/RelWithDebInfo/libglad.a


PostBuild.glfw.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/glfw/src/RelWithDebInfo/libglfw3.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/glfw/src/RelWithDebInfo/libglfw3.a


PostBuild.imgui.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/imgui/RelWithDebInfo/libimgui.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/imgui/RelWithDebInfo/libimgui.a


PostBuild.lexers.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/liblexers.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/liblexers.a


PostBuild.math.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libmath.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libmath.a


PostBuild.predicates.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/predicates/RelWithDebInfo/libpredicates.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/predicates/RelWithDebInfo/libpredicates.a


PostBuild.simd.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libsimd.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libsimd.a


PostBuild.sys.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libsys.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libsys.a


PostBuild.tasking.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libtasking.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libtasking.a


PostBuild.test_bin.RelWithDebInfo:
PostBuild.bcclean.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.embed.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.json.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.igl_cgal.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.CGAL_Core.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.CGAL.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.igl_embree.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.embree.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.simd.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.lexers.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.sys.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.math.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.tasking.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.igl_opengl_glfw_imgui.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.igl_opengl_glfw.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.igl_opengl.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.imgui.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.glfw.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.glad.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.igl_predicates.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.igl.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.igl_common.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.predicates.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.cxxopts.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.nanoflann.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
PostBuild.pybind11.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin
/Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin:\
	/Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/libbcclean.a\
	/usr/local/lib/libboost_thread-mt.a\
	/usr/local/lib/libboost_system-mt.a\
	/usr/local/lib/libboost_chrono-mt.a\
	/usr/local/lib/libboost_date_time-mt.a\
	/usr/local/lib/libboost_atomic-mt.a\
	/usr/local/lib/libgmp.dylib\
	/usr/local/lib/libmpfr.dylib\
	/usr/local/lib/libboost_thread-mt.a\
	/usr/local/lib/libboost_chrono-mt.a\
	/usr/local/lib/libboost_date_time-mt.a\
	/usr/local/lib/libboost_atomic-mt.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libembree3.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libsimd.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/liblexers.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libsys.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libmath.a\
	/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libtasking.a\
	/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/System/Library/Frameworks/OpenGL.framework/OpenGL.tbd\
	/Users/vector_cat/gits/bcClean/Xcode/imgui/RelWithDebInfo/libimgui.a\
	/Users/vector_cat/gits/bcClean/Xcode/glfw/src/RelWithDebInfo/libglfw3.a\
	/Users/vector_cat/gits/bcClean/Xcode/glad/RelWithDebInfo/libglad.a\
	/Users/vector_cat/gits/bcClean/Xcode/predicates/RelWithDebInfo/libpredicates.a\
	/usr/local/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7m.dylib
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/test_bin




# For each target create a dummy ruleso the target does not have to exist
/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/System/Library/Frameworks/OpenGL.framework/OpenGL.tbd:
/Users/vector_cat/gits/bcClean/Xcode/Debug/libbcclean.a:
/Users/vector_cat/gits/bcClean/Xcode/MinSizeRel/libbcclean.a:
/Users/vector_cat/gits/bcClean/Xcode/RelWithDebInfo/libbcclean.a:
/Users/vector_cat/gits/bcClean/Xcode/Release/libbcclean.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libembree3.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/liblexers.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libmath.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libsimd.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libsys.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libtasking.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libembree3.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/liblexers.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libmath.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libsimd.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libsys.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libtasking.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libembree3.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/liblexers.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libmath.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libsimd.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libsys.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libtasking.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libembree3.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/liblexers.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libmath.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libsimd.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libsys.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libtasking.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/bcclean.build/Debug/algorithms.build/Objects-normal/libalgorithms.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/bcclean.build/MinSizeRel/algorithms.build/Objects-normal/libalgorithms.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/bcclean.build/RelWithDebInfo/algorithms.build/Objects-normal/libalgorithms.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/bcclean.build/Release/algorithms.build/Objects-normal/libalgorithms.a:
/Users/vector_cat/gits/bcClean/Xcode/glad/Debug/libglad.a:
/Users/vector_cat/gits/bcClean/Xcode/glad/MinSizeRel/libglad.a:
/Users/vector_cat/gits/bcClean/Xcode/glad/RelWithDebInfo/libglad.a:
/Users/vector_cat/gits/bcClean/Xcode/glad/Release/libglad.a:
/Users/vector_cat/gits/bcClean/Xcode/glfw/src/Debug/libglfw3.a:
/Users/vector_cat/gits/bcClean/Xcode/glfw/src/MinSizeRel/libglfw3.a:
/Users/vector_cat/gits/bcClean/Xcode/glfw/src/RelWithDebInfo/libglfw3.a:
/Users/vector_cat/gits/bcClean/Xcode/glfw/src/Release/libglfw3.a:
/Users/vector_cat/gits/bcClean/Xcode/imgui/Debug/libimgui.a:
/Users/vector_cat/gits/bcClean/Xcode/imgui/MinSizeRel/libimgui.a:
/Users/vector_cat/gits/bcClean/Xcode/imgui/RelWithDebInfo/libimgui.a:
/Users/vector_cat/gits/bcClean/Xcode/imgui/Release/libimgui.a:
/Users/vector_cat/gits/bcClean/Xcode/predicates/Debug/libpredicates.a:
/Users/vector_cat/gits/bcClean/Xcode/predicates/MinSizeRel/libpredicates.a:
/Users/vector_cat/gits/bcClean/Xcode/predicates/RelWithDebInfo/libpredicates.a:
/Users/vector_cat/gits/bcClean/Xcode/predicates/Release/libpredicates.a:
/usr/local/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7m.dylib:
/usr/local/lib/libboost_atomic-mt.a:
/usr/local/lib/libboost_chrono-mt.a:
/usr/local/lib/libboost_date_time-mt.a:
/usr/local/lib/libboost_system-mt.a:
/usr/local/lib/libboost_thread-mt.a:
/usr/local/lib/libgmp.dylib:
/usr/local/lib/libmpfr.dylib:
