include("Draw.jl")

using Gtk
using Printf

win = GtkWindow("NetworkVis")
g = GtkGrid()
netLabel = GtkLabel("Network")
netEntry = GtkEntry()
set_gtk_property!(netEntry, :text, "Karate")
n0Label = GtkLabel("n0")
n0Entry = GtkEntry()
set_gtk_property!(n0Entry, :text, "20")
m0Label = GtkLabel("m0")
m0Entry = GtkEntry()
set_gtk_property!(m0Entry, :text, "20")
selectCb = GtkComboBoxText()
push!(selectCb, "Node")
push!(selectCb, "Node && Edge")
set_gtk_property!(selectCb, :active, 1)
addCb = GtkComboBoxText()
push!(addCb, "KeepOne")
push!(addCb, "KeepAll")
set_gtk_property!(addCb, :active, 0)

retLabel = GtkLabel("*****")
retLabel2 = GtkLabel("*****")
retLabel3 = GtkLabel("*****")

drawButton = GtkButton("Draw")
g[1, 1] = netLabel
g[2, 1] = netEntry
g[1, 2] = n0Label
g[2, 2] = n0Entry
g[1, 3] = m0Label
g[2, 3] = m0Entry
g[1, 4] = selectCb
g[2, 4] = addCb
g[1:2, 5] = drawButton
g[1:2, 6] = retLabel
g[1:2, 7] = retLabel2
g[1:2, 8] = retLabel3

set_gtk_property!(g, :column_homogeneous, true)
set_gtk_property!(g, :column_spacing, 15)
set_gtk_property!(g, :row_spacing, 15)
set_gtk_property!(g, :margin_left, 15)
set_gtk_property!(g, :margin_right, 15)
set_gtk_property!(g, :margin_top, 15)
set_gtk_property!(g, :margin_bottom, 15)
push!(win, g)
showall(win)

signal_connect(drawButton, "clicked") do widget
	s1 = get_gtk_property(netEntry, :text, String)
	s2 = get_gtk_property(n0Entry, :text, String)
	s3 = get_gtk_property(m0Entry, :text, String)
	s4 = Gtk.bytestring( GAccessor.active_text(selectCb) )
	s5 = Gtk.bytestring( GAccessor.active_text(addCb) )
	cntstr = "Network = " * s1 * ", n0 = " * s2 * ", m0 = " * s3
	cntstr2 = "Compress Mode : " * s4 * ", " * s5
	GAccessor.text(retLabel, cntstr)
	GAccessor.text(retLabel2, cntstr2)
	GAccessor.text(retLabel3, "Running ... ...")
	ags1 = s1
	ags2 = s2 * "," * s3 * "," * ((s4 == "Node") ? "n," : "ne,") * s5
	T = time()
	DrawNetwork(ags1, ags2)
	T = time() - T
	GAccessor.text(retLabel3, "Running Time : " * @sprintf("%.2f s", T))
	run(Cmd(["epstopdf", "../pictures/"*ags1*".eps"]))
	run(Cmd(["open", "../pictures/"*ags1*".pdf"]))

	#win2 = GtkWindow("Final Drawing")
	#bx = GtkBox(:v)
	#ii = GtkImage("picture/" * s1 * ".eps")
	#push!(bx, ii)
	#push!(win2, bx)
	#showall(win2)
	#signal_connect(win2, :destroy) do widget
	#    Gtk.gtk_quit()
	#end
end

signal_connect(win, :destroy) do widget
    Gtk.gtk_quit()
	exit()
end

Gtk.gtk_main()
