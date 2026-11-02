{pkgs ? import <nixpkgs> {}}:
pkgs.mkShell {
  buildInputs = with pkgs; [
    # Build systems and tools
    cmake
    pkg-config
    gnumake
    gcc

    # OpenGL Core & Utilities
    libGL
    libGLU
    glew
    glfw

    # X11 Core & Extensions (Fixes the missing -lXxf86vm, -lXrandr, and -lXi errors)
    libX11
    libX11.dev
    xorgproto
    libXxf86vm
    libXrandr
    libXi
    libXinerama # Adding this as well, as GLFW usually requires it next
    libXcursor # Adding this as well for mouse cursor management
  ];

  shellHook = ''
    export LD_LIBRARY_PATH=/run/opengl-driver/lib:/run/opengl-driver-32/lib:$LD_LIBRARY_PATH
  '';
}
