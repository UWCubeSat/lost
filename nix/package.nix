{ stdenv, pkgs, lib }:

stdenv.mkDerivation rec {
  pname = "lost";
  version = "0.1.0";

  src = ./..;

  nativeBuildInputs = with pkgs; [ gcc gnumake groff unixtools.xxd eigen ];

  buildInputs = with pkgs; [ cairo ];

  dontConfigure = true;

  buildPhase = ''
    make release
  '';

  installPhase = ''
    mkdir -p $out/bin
    mv lost $out/bin/
  '';

  meta = with lib; {
    description = "Open-source Star Tracker";
    homepage = "https://github.com/uwcubesat/lost";
    license = licenses.mit;
  };
}
