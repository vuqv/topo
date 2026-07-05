# view_trunc.tcl — visualize a truncated CG ribosome and draw the coordinate axes
# to confirm the exit tunnel is aligned on the X-axis (tunnel line ≈ (x, 0, 0),
# PTC at the origin, exit toward +x).
#
# Run:
#   vmd -e assets/csp/prepare_ribosome/helpers/view_trunc.tcl                 # uses $pdbfile below
#   vmd yourfile.pdb -e assets/csp/prepare_ribosome/helpers/view_trunc.tcl    # or load first, script skips its own load
#   # or inside VMD:  source assets/csp/prepare_ribosome/helpers/view_trunc.tcl
#
# What you should see if orientation is correct:
#   * the RED X-axis threads straight down the empty tunnel channel,
#   * the beads near x≈0 (the PTC / tRNA CCA ends) sit tight around the axis,
#   * the structure opens up toward +x (the cytosolic exit).
# ---------------------------------------------------------------------------

# ---- user-editable parameters ---------------------------------------------
# Path relative to the bundle root (assets/csp/prepare_ribosome/); any shipped
# structure works. Run e.g.:  vmd -e helpers/view_trunc.tcl
set pdbfile "structures/human/8g61_60S_model_cg_trunc.pdb"

# Axis extents (Å). X spans the whole structure; Y/Z are a short reference triad.
set x_min  -25.0
set x_max  120.0
set yz_len  40.0
# ---------------------------------------------------------------------------

# ---- load structure (skip if VMD already has one) -------------------------
if {[molinfo num] == 0} {
    if {![file exists $pdbfile]} {
        error "pdbfile not found: $pdbfile (edit the path at the top of the script)"
    }
    mol new $pdbfile type pdb waitfor all
}
set molid [molinfo top]
mol delrep 0 $molid

# ---- representations -------------------------------------------------------
# All beads, small VDW, colored by segID (VMD maps PDB segID → SegName).
mol representation VDW 0.6 12.0
mol color SegName
mol selection {all}
mol material Opaque
mol addrep $molid

# Highlight the two tRNAs (near the PTC) so you can see they sit at x≈0.
mol representation VDW 1.2 16.0
mol color ColorID 3   ;# orange
mol selection {segname PtR}
mol addrep $molid

mol representation VDW 1.2 16.0
mol color ColorID 11  ;# purple
mol selection {segname AtR}
mol addrep $molid

# ---- draw the axes as arrows (X red, Y green, Z blue) ----------------------
proc draw_arrow {mol start end colorid {rad 0.8}} {
    graphics $mol color $colorid
    set shaft [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
    graphics $mol cylinder $start $shaft radius $rad resolution 24 filled yes
    graphics $mol cone $shaft $end radius [expr {$rad * 2.5}] resolution 24
}

# clear any previous drawings on this molecule
graphics $molid delete all

# X-axis (red = ColorID 1): the exit-tunnel axis. This is the one to check.
draw_arrow $molid [list $x_min 0 0] [list $x_max 0 0] 1 0.8
# Y-axis (green = 7) and Z-axis (blue = 0): short reference triad at the origin.
draw_arrow $molid [list 0 0 0] [list 0 $yz_len 0] 7 0.6
draw_arrow $molid [list 0 0 0] [list 0 0 $yz_len] 0 0.6

# Origin marker = PTC (yellow sphere).
graphics $molid color 4
graphics $molid sphere {0 0 0} radius 2.5 resolution 24

# Axis labels.
graphics $molid color 1
graphics $molid text [list [expr {$x_max + 3}] 0 0] "X (tunnel/exit)" size 1.2 thickness 2
graphics $molid color 7
graphics $molid text [list 0 [expr {$yz_len + 3}] 0] "Y" size 1.2 thickness 2
graphics $molid color 0
graphics $molid text [list 0 0 [expr {$yz_len + 3}]] "Z" size 1.2 thickness 2

# ---- display settings ------------------------------------------------------
color Display Background white
display projection Orthographic
display depthcue off
display resetview
axes location off   ;# hide VMD's corner axes; we drew our own at the origin

puts "Loaded [molinfo $molid get name]"
puts "Drew X (red, tunnel), Y (green), Z (blue); origin sphere = PTC."
puts "Tunnel is correct if the red X-axis runs down the channel and PTC/tRNA beads hug x≈0."
