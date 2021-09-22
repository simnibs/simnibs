from simnibs import SIMNIBSDIR
import os

simnibs_logo = os.path.join(SIMNIBSDIR, '_internal_resources', 'icons', 'simnibs', 'gui_icon.png')

html_source = f'''
<!DOCTYPE html>

<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Charm summary</title>
</head>

<style>
    #content {
        position: relative;
    }
    #content img {
        position: absolute;
        top: 0px;
        right: 0px;
    }
</style>

<div id="content">
    <img src="%s" class="ribbon"/>
</div>

<body>
<h1> Charm run summary </h1>
<br><a href="reg_viewer.html">T1-T2 registration</a>
<a href="final_viewer.html">Final segmentation</a>

<h1> Quality check </h1>
<a href="reg_viewer.html">T1-T2 registration</a>
<br>
<br> <a href="final_viewer.html">Final segmentation</a>

<h1> References </h1>
<a href="https://www.sciencedirect.com/science/article/pii/S1053811920305309"> Puonti, O., et al., "Accurate and robust whole-head segmentation from magnetic resonance images for individualized head modeling", NeuroImage, 219, 2020 </a>
<br>
<br> <a href="https://ieeexplore.ieee.org/document/7318340"> Thielscher, A., et al., "Field modeling for transcranial magnetic stimulation: a useful tool to understand the physiological effects of TMS?", IEEE EMBS 2015, Milano, Italy </a>
</body>
</html>
'''

def write_template():

text

