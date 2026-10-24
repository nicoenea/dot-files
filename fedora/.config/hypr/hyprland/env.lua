------------------
---- GPU ENVS ----
------------------
-- AMD RX 9070 XT (amdgpu/RADV). Do NOT reintroduce __GLX_VENDOR_LIBRARY_NAME or
-- NVD_BACKEND — those are NVIDIA-only and break GLX/VA-API once the driver is gone.
local home = os.getenv("HOME")

hl.env("LIBVA_DRIVER_NAME", "radeonsi")
hl.env("ELECTRON_OZONE_PLATFORM_HINT", "auto")
hl.env("GDK_BACKEND", "wayland,x11,*")
hl.env("QT_QPA_PLATFORM", "wayland;xcb")
hl.env("SDL_VIDEODRIVER", "wayland")
hl.env("CLUTTER_BACKEND", "wayland")
hl.env("XDG_CURRENT_DESKTOP", "Hyprland")
hl.env("XDG_SESSION_TYPE", "wayland")
hl.env("XDG_SESSION_DESKTOP", "Hyprland")
-- caelestia-shell/m3shapes are built into /usr/local by update-caelestia;
-- Qt and the shell only look under /usr by default.
hl.env("QML_IMPORT_PATH", "/usr/local/lib64/qt6/qml")
hl.env("CAELESTIA_LIB_DIR", "/usr/local/lib64/caelestia")
hl.env("QT_QPA_PLATFORMTHEME", "gtk3")
hl.env("PATH", home .. "/.local/bin:/usr/local/bin:/usr/bin:/bin")
