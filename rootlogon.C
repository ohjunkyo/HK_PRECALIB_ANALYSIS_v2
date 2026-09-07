// Auto-loaded by ROOT whenever a session starts in this directory (root -l,
// root -l -b -q, etc. all read rootlogon.C from cwd unless -n is passed).
//
// ACLiC (.x Macro.C+) used to compile its .so/.d/.pcm output right next to
// the source file, littering this directory with per-macro build artifacts
// that don't belong in the analysis repo. This redirects that build cache to
// a dotdir here instead -- still local (fast rebuilds, no cross-user
// clashes), just out of the way and easy to wipe (`rm -rf .aclic_cache`).
{
    gSystem->mkdir(".aclic_cache", kTRUE);
    gSystem->SetBuildDir(".aclic_cache", kTRUE);
}
