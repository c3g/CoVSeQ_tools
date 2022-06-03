# Here are the release instructions

# Version bump the value. Remove '-beta'
vim VERSION

# Update mugqic_tools module version
vim ../genpipes/resources/modules/covseq_tools.sh
(VERSION=1.2.1)

# Tag the branch and push the tag. You'll need to have a gpg signature for this. Extra precaution
git tag 1.2.1 -m 'Release 1.2.1'
git push -u origin --tags

# Recreate the CHANGELOG.md
bash ~/repo/dump_ChangeLog.sh > CHANGELOG.md
git commit -a -m "Version bump to 1.2.1"
git push

# Create a release tarball archive
git archive --format=tar --prefix=covseq_tools-1.2.1/ <latest_commmit> | gzip > ~/covseq_tools-1.2.1.tar.gz

# Upload this archive in
https://bitbucket.org/mugqic/covseq_tools/downloads

# Version bump the value. Until the next release, add '-beta' e.g. 1.2.2-beta
vim VERSION
git commit -m "Version bump to 1.2.2-beta" VERSION
git push

# Deploy covseq_tools-<VERSION> as a module on all clusters
