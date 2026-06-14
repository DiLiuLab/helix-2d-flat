# -*- mode: python ; coding: utf-8 -*-

import os
import sys


project_root = os.path.abspath(SPECPATH)
png_icon = os.path.join(project_root, "assets", "helix_2d_flat_icon.png")
if sys.platform == "darwin":
    executable_icon = os.path.join(project_root, "assets", "helix_2d_flat_icon.icns")
elif sys.platform == "win32":
    executable_icon = os.path.join(project_root, "assets", "helix_2d_flat_icon.ico")
else:
    executable_icon = None

a = Analysis(
    ["helix_2D_flatV3.py"],
    pathex=[project_root],
    binaries=[],
    datas=[(png_icon, "assets")],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

if sys.platform == "darwin":
    exe = EXE(
        pyz,
        a.scripts,
        [],
        exclude_binaries=True,
        name="Helix2DFlat",
        debug=False,
        bootloader_ignore_signals=False,
        strip=False,
        upx=True,
        console=False,
        argv_emulation=False,
        target_arch=None,
        codesign_identity=None,
        entitlements_file=None,
        icon=executable_icon,
    )
    collected = COLLECT(
        exe,
        a.binaries,
        a.datas,
        strip=False,
        upx=True,
        upx_exclude=[],
        name="Helix2DFlat",
    )
    app = BUNDLE(
        collected,
        name="Helix 2D Flat.app",
        icon=executable_icon,
        bundle_identifier="org.diliulab.helix-2d-flat",
        info_plist={
            "CFBundleDisplayName": "Helix 2D Flat",
            "NSHighResolutionCapable": True,
        },
    )
else:
    exe = EXE(
        pyz,
        a.scripts,
        a.binaries,
        a.datas,
        [],
        name="Helix2DFlat",
        debug=False,
        bootloader_ignore_signals=False,
        strip=False,
        upx=True,
        upx_exclude=[],
        runtime_tmpdir=None,
        console=False,
        disable_windowed_traceback=False,
        icon=executable_icon,
    )
