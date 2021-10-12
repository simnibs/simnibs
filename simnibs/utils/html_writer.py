from simnibs import SIMNIBSDIR
import os
import json
from simnibs.utils.settings_reader import read_ini


'''
This module writes the final html after a charm run is done.
The basic .html template is defined below, which is
pieced together from different text snippets and links
depending on what charm was run on.
'''

simnibs_logo = os.path.join(SIMNIBSDIR,
                            '_internal_resources',
                            'icons', 'simnibs',
                            'gui_icon.png')

html_source = """
<!DOCTYPE html>

<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Charm summary</title>
  <style>
  div {{
    max-width:1000px;
    min-width:100px;
    align:center;
    font-family:sans-serif;
  }}
</style>
</head>

<div id="content">
    <a href="https://simnibs.github.io/simnibs/build/html/index.html">
    <img src="{simnibs_fig}" width="128" height="128" class="ribbon"/> </a>
</div>


<body>
<h2> Charm run summary: </h2>
<div>
{scan_text} {settings_text}
Please check the quality of the inital registration
between the input scan and the atlas by clicking the
affine viewer link below,
and the quality of the final
segmentation by clicking the final segmentation viewer
link below.
The registration between the T1- and
T2-weighted scans can be inspected by clicking the
registration viewer link.
Finally, if you use charm and/or
SimNIBS in your study, please include citations to the papers
listed under the "References"-section
</div>

<h2> Viewer links for quality checking </h2>
<div>
{t1_t2_reg}
<br>
<br>{affine_reg}
<br>
<br>{final_seg}
</div>

<h2> References </h2>
<div>
<a href="https://www.sciencedirect.com/science/article/pii/S1053811920305309">
Puonti, O., et al., "Accurate and robust whole-head segmentation from magnetic
resonance images for individualized head modeling", NeuroImage, 219, 2020 </a>
<br>
<br> <a href="https://ieeexplore.ieee.org/document/7318340">
Thielscher, A., et al., "Field modeling for transcranial magnetic stimulation:
a useful tool to understand the physiological effects of TMS?",
IEEE EMBS 2015, Milano, Italy </a>
</div>
</body>
</html>
"""


def _explode_nested_list(x, output=None):
    """Explode a list.

    Flattens an arbitrarily nested list.
    Returns a flattened list.

    Parameters
    ----------------------------------
    x: list
    output: flattened list
    """
    if output is None:
        output = []

    for i in x:
        if isinstance(i, list):
            _explode_nested_list(i, output)
        else:
            output.append(i)

    return output


def _is_primitive(x):
    """Check if the input is primitive.

    Check if input is a "primitive" type.
    Returns True if x is int, float, str, chr or boolean,
    otherwise False.

    Parameters
    ----------------------------------
    x: input variable to check
    """
    try:
        iter(x)
        return False
    except TypeError:
        primitive_types = (int, float, str, bool, chr)
        return isinstance(x, primitive_types)


def _walk_through_dict_recursively(dict1, dict2, diff_dict=None):
    """Compare one dict to another.

    Takes in a dictionary and recursively loops through it,
    and compares the elements (not very specifically).
    Returns those dict1 elements which are not equal to dict2 elements.

    Parameters
    ------------------------------------
    dict1: dictionary
    dict2: dictionary
    diff_dict: dictionary of differing elements
    """
    if diff_dict is None:
        diff_dict = {}

    for key1, key2 in zip(dict1, dict2):

        if isinstance(dict1[key1], dict) and isinstance(dict2[key2], dict):

            if dict1[key1] == dict2[key2]:
                continue

            else:
                diff_dict[key1] = {}
                _walk_through_dict_recursively(dict1[key1],
                                               dict2[key2],
                                               diff_dict[key1])

        elif ((not isinstance(dict1[key1], dict) and
               isinstance(dict2[key2], dict)) or
              (isinstance(dict1[key1], dict) and
               not isinstance(dict2[key2], dict))):
            diff_dict[key1] = dict1[key1]
            continue

        elif isinstance(dict1[key1], list) and isinstance(dict2[key2], list):

            flat_list1 = _explode_nested_list(dict1[key1])
            flat_list2 = _explode_nested_list(dict2[key2])

            if set(flat_list1) == set(flat_list2):
                continue
            else:
                diff_dict[key1] = dict1[key1]

        elif ((isinstance(dict1[key1], list) and
               not isinstance(dict2[key2], list)) or
              (not isinstance(dict1[key1], list) and
               isinstance(dict2[key2], list))):

            diff_dict[key1] = dict1[key1]

        elif _is_primitive(dict1[key1]) and _is_primitive(dict2[key2]):

            if dict1[key1] == dict2[key2]:
                continue
            else:
                diff_dict[key1] = dict1[key1]

        else:
            continue

    return diff_dict


def _get_settings_string(sub_ini, charm_ini):
    """Return difference between two dicts as string.

    Compare to ini files, find fields that differ,
    and return the differing fields in a string.

    Parameters
    -------------------------------
    sub_ini: str

    Path to subject ini file.

    charm_ini: str

    Path to the charm.ini file defining the standard settings.
    """
    sub_settings = read_ini(sub_ini)
    charm_settings = read_ini(charm_ini)
    diff_dict = _walk_through_dict_recursively(sub_settings, charm_settings)

    if not bool(diff_dict):
        settings_text = 'using the standard settings listed in the charm.ini file.<br><br>'
    else:
        settings_text = 'using custom settings for the following fields <br><br>'
        settings_text += json.dumps(diff_dict,
                                    indent=4)

        settings_text += ' <br><br>'
        settings_text = settings_text.replace('\n', '\n<br>')
    return settings_text


def write_template(sub_files):
    """
    Parse and write the .html template to disk.

    Parameters
    -----------------------
    sub_files: SubjectFiles object

    Instantiation of the SubjectFiles class
    for the given subject folder.
    """
    t1_t2_reg_viewer = '<a href="' + sub_files.t1_t2_reg_viewer + '">T1-T2 registration viewer</a>'

    aff_reg_viewer = '<a href="' + sub_files.affine_reg_viewer + '">Affine registration viewer</a>'

    final_seg_viewer = '<a href="' + sub_files.final_seg_viewer + '">Final segmentation viewer</a>'

    t1_text = 'Charm was run on a T1-weighted scan'
    t1_t2_text = 'Charm was run on a combination of T1- and T2-weighted scans'

    parse_dict = {}
    if not os.path.exists(sub_files.T2_reg):
        parse_dict['scan_text'] = t1_text
    else:
        parse_dict['scan_text'] = t1_t2_text

    parse_dict['simnibs_fig'] = simnibs_logo

    if not os.path.exists(sub_files.t1_t2_reg_viewer):
        parse_dict["t1_t2_reg"] = ''
    else:
        parse_dict["t1_t2_reg"] = t1_t2_reg_viewer

    if not os.path.exists(sub_files.affine_reg_viewer):
        parse_dict["affine_reg"] = ''
    else:
        parse_dict["affine_reg"] = aff_reg_viewer

    if not os.path.exists(sub_files.final_seg_viewer):
        parse_dict["final_seg"] = ''
    else:
        parse_dict["final_seg"] = final_seg_viewer

    settings_dif = _get_settings_string(sub_files.settings,
                                        os.path.join(SIMNIBSDIR, 'charm.ini'))
    parse_dict["settings_text"] = settings_dif

    html_parsed = html_source.format(**parse_dict)
    with open(sub_files.summary_report, "w") as f:
        f.write(html_parsed)
