gfx create region mesh
gfx read region reacdiff.xml region mesh
gfx create window
gfx define faces egroup mesh
gfx edit scene

gfx modify g_element /mesh/ general clear;
gfx modify g_element /mesh/ lines coordinate reacdiff.geometric tessellation default LOCAL select_on material default selected_material default_selected;
gfx modify g_element /mesh/ element_points coordinate reacdiff.geometric discretization "1*1*1" tessellation NONE LOCAL glyph point general size "1*1*1" centre 0,0,0 font default use_elements cell_density density reacdiff.dependent select_on material default data reacdiff.dependent spectrum default selected_material default_selected;
gfx edit spectrum
gfx list node region /mesh 88
