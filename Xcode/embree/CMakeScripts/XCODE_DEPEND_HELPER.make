# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.algorithms.Debug:
PostBuild.embree.Debug:
PostBuild.algorithms.Debug: /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libembree3.a
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libembree3.a:\
	/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/embree3.build/Debug/algorithms.build/Objects-normal/libalgorithms.a
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libembree3.a


PostBuild.lexers.Debug:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/liblexers.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/liblexers.a


PostBuild.math.Debug:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libmath.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libmath.a


PostBuild.simd.Debug:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libsimd.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libsimd.a


PostBuild.sys.Debug:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libsys.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libsys.a


PostBuild.tasking.Debug:
/Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libtasking.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Debug/libtasking.a


PostBuild.algorithms.Release:
PostBuild.embree.Release:
PostBuild.algorithms.Release: /Users/vector_cat/gits/bcClean/Xcode/embree/Release/libembree3.a
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libembree3.a:\
	/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/embree3.build/Release/algorithms.build/Objects-normal/libalgorithms.a
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Release/libembree3.a


PostBuild.lexers.Release:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/liblexers.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Release/liblexers.a


PostBuild.math.Release:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libmath.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Release/libmath.a


PostBuild.simd.Release:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libsimd.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Release/libsimd.a


PostBuild.sys.Release:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libsys.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Release/libsys.a


PostBuild.tasking.Release:
/Users/vector_cat/gits/bcClean/Xcode/embree/Release/libtasking.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/Release/libtasking.a


PostBuild.algorithms.MinSizeRel:
PostBuild.embree.MinSizeRel:
PostBuild.algorithms.MinSizeRel: /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libembree3.a
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libembree3.a:\
	/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/embree3.build/MinSizeRel/algorithms.build/Objects-normal/libalgorithms.a
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libembree3.a


PostBuild.lexers.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/liblexers.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/liblexers.a


PostBuild.math.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libmath.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libmath.a


PostBuild.simd.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libsimd.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libsimd.a


PostBuild.sys.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libsys.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libsys.a


PostBuild.tasking.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libtasking.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/MinSizeRel/libtasking.a


PostBuild.algorithms.RelWithDebInfo:
PostBuild.embree.RelWithDebInfo:
PostBuild.algorithms.RelWithDebInfo: /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libembree3.a
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libembree3.a:\
	/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/embree3.build/RelWithDebInfo/algorithms.build/Objects-normal/libalgorithms.a
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libembree3.a


PostBuild.lexers.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/liblexers.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/liblexers.a


PostBuild.math.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libmath.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libmath.a


PostBuild.simd.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libsimd.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libsimd.a


PostBuild.sys.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libsys.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libsys.a


PostBuild.tasking.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libtasking.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/embree/RelWithDebInfo/libtasking.a




# For each target create a dummy ruleso the target does not have to exist
/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/embree3.build/Debug/algorithms.build/Objects-normal/libalgorithms.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/embree3.build/MinSizeRel/algorithms.build/Objects-normal/libalgorithms.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/embree3.build/RelWithDebInfo/algorithms.build/Objects-normal/libalgorithms.a:
/Users/vector_cat/gits/bcClean/Xcode/embree/common/algorithms/embree3.build/Release/algorithms.build/Objects-normal/libalgorithms.a:
