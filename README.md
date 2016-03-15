ProbeDeformer plugins for Maya
/**
 * @brief Probe Deformer plugin for Maya
 * @section LICENSE The MIT License
 * @section requirements:  Eigen 3:  http://eigen.tuxfamily.org/
 * @section Autodesk Maya: http://www.autodesk.com/products/autodesk-maya/overview
 * @section (included) AffineLib: https://github.com/shizuo-kaji/AffineLib
 * @version 0.10
 * @date  3/Nov/2013
 * @author Shizuo KAJI
 */

There are two versions of deformers:
1. "probeDeformer" is a simple vertex based deformer which ignores mesh structure.
It is applicable to any mesh and particles.
2. "probeDeformerARAP" is a mesh based deformer which try to preserve geometry (ARAP).
It is applicable only to "clean" meshes; it is recommended first "cleanup" the target mesh by
"Cleanup" => "Remove Zero edges, faces" in Maya's mesh menu.

For the detail of the algorithm, refer to the papers
* "Probe-Type Deformers" by S. Kaji and G. Liu
in Mathematical Progress in Expressive Image Synthesis II, Volume 18 of the series Mathematics for Industry pp 63-77
http://link.springer.com/chapter/10.1007/978-4-431-55483-7_6?no-access=true
* "A concise parametrisation of affine transformation" by S. Kaji and H. Ochiai.
http://arxiv.org/abs/1507.05290

How to compile:
* For Mac users, look at the included Xcode project file ( or Makefile )
* For Windows users, look at the included Visual Studio project file.
Please refer to Autodesk's web page for details.

How to use:
1. Place the plugin files in "MAYA_PLUG_IN_PATH"
2. Place the UI python script files in "MAYA_SCRIPT_PATH"
3. Open Script editor in Maya and type in the following Python command:

    import ui_probeDeformer as ui
    ui.UI_ProbeDeformer()


To visualise vertex color, go to "Display" => "Polygon" => "Custom Polygon Display"
and tick "color" and select "emission."

LIMITATION:
The ARAP version works only on "clean" meshes.
First apply "Cleanup" from "Mesh" menu
to remove zero area faces and zero length edges.

