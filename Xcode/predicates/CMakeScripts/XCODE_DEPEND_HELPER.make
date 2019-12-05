# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.predicates.Debug:
/Users/vector_cat/gits/bcClean/Xcode/predicates/Debug/libpredicates.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/predicates/Debug/libpredicates.a


PostBuild.predicates.Release:
/Users/vector_cat/gits/bcClean/Xcode/predicates/Release/libpredicates.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/predicates/Release/libpredicates.a


PostBuild.predicates.MinSizeRel:
/Users/vector_cat/gits/bcClean/Xcode/predicates/MinSizeRel/libpredicates.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/predicates/MinSizeRel/libpredicates.a


PostBuild.predicates.RelWithDebInfo:
/Users/vector_cat/gits/bcClean/Xcode/predicates/RelWithDebInfo/libpredicates.a:
	/bin/rm -f /Users/vector_cat/gits/bcClean/Xcode/predicates/RelWithDebInfo/libpredicates.a




# For each target create a dummy ruleso the target does not have to exist
