#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR/../..

IFS=$'\n'
filter='opt/local/include/\|opt/local/lib/\|opt/local/share/man/\|Applications/Embree3'
FILES+=($(pkgutil --files com.intel.embree-3.5.2.examples | grep -e $filter | tail -r))
FILES+=($(pkgutil --files com.intel.embree-3.5.2.lib      | grep -e $filter | tail -r))
FILES+=($(pkgutil --files com.intel.embree-3.5.2.devel    | grep -e $filter | tail -r))
unset IFS

# exit if no files found
if [ ${#FILES[@]} -eq 0 ]; then
  printf "Embree 3.5.2 not installed!\n"
  exit
fi

# first print all files that would get removed
echo Uninstalling Embree 3.5.2 will remove the following files:
PWD=`pwd`
if [ "$PWD" != "/" ]; then
    PWD=$PWD/
fi
for f in "${FILES[@]}"; do
    printf "  %s%s\n" $PWD "$f"
done

echo "Do you wish to uninstall Embree 3.5.2 by removing these files?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) break;;
        No ) exit;;
    esac
done

# now remove files
echo Uninstalling Embree 3.5.2 ...
for f in "${FILES[@]}"; do
    sudo /bin/rm -vd "$f"
done

sudo /usr/sbin/pkgutil --forget com.intel.embree-3.5.2.examples
sudo /usr/sbin/pkgutil --forget com.intel.embree-3.5.2.devel
sudo /usr/sbin/pkgutil --forget com.intel.embree-3.5.2.lib
