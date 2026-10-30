{
  description = "Open-source Star Tracker";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/release-24.05";
    flake-utils.url = "github:numtide/flake-utils";

    flake-compat = {
      url = "github:edolstra/flake-compat";
      flake = false;
    };
  };

  outputs = inputs@{ self, nixpkgs, flake-utils, ... }:
    let
      overlays = [
        (final: prev: rec { lost = prev.callPackage ./nix/package.nix { }; })
      ];
    in flake-utils.lib.eachDefaultSystem (system:
      let pkgs = import nixpkgs { inherit overlays system; };
      in {
        packages.default = pkgs.lost;
        packages.lost = pkgs.lost;

        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            gcc
            clang-tools
            gnumake
            bear
            gdb
            valgrind
            man
            groff
            imagemagick
            cpplint
            doxygen
            graphviz
            unixtools.xxd
            eigen
            cairo
          ];
        };

        devShell = self.devShells.${system}.default;
      });
}
