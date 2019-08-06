File = "snap.py"
filetag(v) = Float64(Int(round(v*100000))) / 100000
filename(i) = begin
	file = "snap"
	if i < 10
		file = file * "000$(i)"
	elseif i < 100
		file = file * "00$(i)"
	end
end
lambdas = collect(linearspace(0.0, 1.0, 49))
for i in 1:10
	push!(lambdas, lambdas[end])
end
for (i, lambda) in enumerate(lambdas[2:end])

open(File, "w") do f
	print(f, """
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
bertoldi_compressionboundaryvtk = LegacyVTKReader(FileNames=['bertoldi_compression-boundary.vtk'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [966, 966]

# show data in view
bertoldi_compressionboundaryvtkDisplay = Show(bertoldi_compressionboundaryvtk, renderView1)

# trace defaults for the display properties.
bertoldi_compressionboundaryvtkDisplay.Representation = 'Surface'
bertoldi_compressionboundaryvtkDisplay.ColorArrayName = [None, '']
bertoldi_compressionboundaryvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
bertoldi_compressionboundaryvtkDisplay.SelectOrientationVectors = 'None'
bertoldi_compressionboundaryvtkDisplay.ScaleFactor = 0.0101258232648
bertoldi_compressionboundaryvtkDisplay.SelectScaleArray = 'None'
bertoldi_compressionboundaryvtkDisplay.GlyphType = 'Arrow'
bertoldi_compressionboundaryvtkDisplay.GlyphTableIndexArray = 'None'
bertoldi_compressionboundaryvtkDisplay.GaussianRadius = 0.00050629116324
bertoldi_compressionboundaryvtkDisplay.SetScaleArray = [None, '']
bertoldi_compressionboundaryvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
bertoldi_compressionboundaryvtkDisplay.OpacityArray = [None, '']
bertoldi_compressionboundaryvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
bertoldi_compressionboundaryvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
bertoldi_compressionboundaryvtkDisplay.SelectionCellLabelFontFile = ''
bertoldi_compressionboundaryvtkDisplay.SelectionPointLabelFontFile = ''
bertoldi_compressionboundaryvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
bertoldi_compressionboundaryvtkDisplay.ScalarOpacityUnitDistance = 0.011328756643810481

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
bertoldi_compressionboundaryvtkDisplay.DataAxesGrid.XTitleFontFile = ''
bertoldi_compressionboundaryvtkDisplay.DataAxesGrid.YTitleFontFile = ''
bertoldi_compressionboundaryvtkDisplay.DataAxesGrid.ZTitleFontFile = ''
bertoldi_compressionboundaryvtkDisplay.DataAxesGrid.XLabelFontFile = ''
bertoldi_compressionboundaryvtkDisplay.DataAxesGrid.YLabelFontFile = ''
bertoldi_compressionboundaryvtkDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
bertoldi_compressionboundaryvtkDisplay.PolarAxes.PolarAxisTitleFontFile = ''
bertoldi_compressionboundaryvtkDisplay.PolarAxes.PolarAxisLabelFontFile = ''
bertoldi_compressionboundaryvtkDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
bertoldi_compressionboundaryvtkDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.CameraPosition = [0.0005937481759999974, 0.0504437499, 10000.0]
renderView1.CameraFocalPoint = [0.0005937481759999974, 0.0504437499, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Legacy VTK Reader'
snapshotvtk = LegacyVTKReader(FileNames=['bertoldi_compression-$(filetag(lambda)).vtk'])

# show data in view
snapshotvtkDisplay = Show(snapshotvtk, renderView1)

# trace defaults for the display properties.
snapshotvtkDisplay.Representation = 'Surface'
snapshotvtkDisplay.ColorArrayName = [None, '']
snapshotvtkDisplay.OSPRayScaleArray = 'u'
snapshotvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
snapshotvtkDisplay.SelectOrientationVectors = 'u'
snapshotvtkDisplay.ScaleFactor = 0.0101258232648
snapshotvtkDisplay.SelectScaleArray = 'None'
snapshotvtkDisplay.GlyphType = 'Arrow'
snapshotvtkDisplay.GlyphTableIndexArray = 'None'
snapshotvtkDisplay.GaussianRadius = 0.00050629116324
snapshotvtkDisplay.SetScaleArray = ['POINTS', 'u']
snapshotvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
snapshotvtkDisplay.OpacityArray = ['POINTS', 'u']
snapshotvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
snapshotvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
snapshotvtkDisplay.SelectionCellLabelFontFile = ''
snapshotvtkDisplay.SelectionPointLabelFontFile = ''
snapshotvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
snapshotvtkDisplay.ScalarOpacityUnitDistance = 0.009992410134288721

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
snapshotvtkDisplay.DataAxesGrid.XTitleFontFile = ''
snapshotvtkDisplay.DataAxesGrid.YTitleFontFile = ''
snapshotvtkDisplay.DataAxesGrid.ZTitleFontFile = ''
snapshotvtkDisplay.DataAxesGrid.XLabelFontFile = ''
snapshotvtkDisplay.DataAxesGrid.YLabelFontFile = ''
snapshotvtkDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
snapshotvtkDisplay.PolarAxes.PolarAxisTitleFontFile = ''
snapshotvtkDisplay.PolarAxes.PolarAxisLabelFontFile = ''
snapshotvtkDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
snapshotvtkDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Warp By Vector'
warpByVector1 = WarpByVector(Input=snapshotvtk)
warpByVector1.Vectors = ['POINTS', 'u']

# show data in view
warpByVector1Display = Show(warpByVector1, renderView1)

# trace defaults for the display properties.
warpByVector1Display.Representation = 'Surface'
warpByVector1Display.ColorArrayName = [None, '']
warpByVector1Display.OSPRayScaleArray = 'u'
warpByVector1Display.OSPRayScaleFunction = 'PiecewiseFunction'
warpByVector1Display.SelectOrientationVectors = 'u'
warpByVector1Display.ScaleFactor = 0.010142636760398721
warpByVector1Display.SelectScaleArray = 'None'
warpByVector1Display.GlyphType = 'Arrow'
warpByVector1Display.GlyphTableIndexArray = 'None'
warpByVector1Display.GaussianRadius = 0.0005071318380199361
warpByVector1Display.SetScaleArray = ['POINTS', 'u']
warpByVector1Display.ScaleTransferFunction = 'PiecewiseFunction'
warpByVector1Display.OpacityArray = ['POINTS', 'u']
warpByVector1Display.OpacityTransferFunction = 'PiecewiseFunction'
warpByVector1Display.DataAxesGrid = 'GridAxesRepresentation'
warpByVector1Display.SelectionCellLabelFontFile = ''
warpByVector1Display.SelectionPointLabelFontFile = ''
warpByVector1Display.PolarAxes = 'PolarAxesRepresentation'
warpByVector1Display.ScalarOpacityUnitDistance = 0.009938686024079718

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
warpByVector1Display.DataAxesGrid.XTitleFontFile = ''
warpByVector1Display.DataAxesGrid.YTitleFontFile = ''
warpByVector1Display.DataAxesGrid.ZTitleFontFile = ''
warpByVector1Display.DataAxesGrid.XLabelFontFile = ''
warpByVector1Display.DataAxesGrid.YLabelFontFile = ''
warpByVector1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
warpByVector1Display.PolarAxes.PolarAxisTitleFontFile = ''
warpByVector1Display.PolarAxes.PolarAxisLabelFontFile = ''
warpByVector1Display.PolarAxes.LastRadialAxisTextFontFile = ''
warpByVector1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(snapshotvtk, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# change representation type
warpByVector1Display.SetRepresentationType('Surface With Edges')

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0012350857876053747, 0.050203248295647984, 0.27629011483940236]
renderView1.CameraFocalPoint = [0.0012350857876053747, 0.050203248295647984, 0.0]
renderView1.CameraParallelScale = 0.05252248397463676

RenderAllViews()
# save screenshot
SaveScreenshot('$(filename(i)).png', renderView1, ImageResolution=[966, 966])

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0012350857876053747, 0.050203248295647984, 0.27629011483940236]
renderView1.CameraFocalPoint = [0.0012350857876053747, 0.050203248295647984, 0.0]
renderView1.CameraParallelScale = 0.05252248397463676

#### uncomment the following to render all views

# alternatively, if you want to write images, you can use SaveScreenshot(...).
quit()
""")
end

run(`"pvpython" $File`)

end

run(`"C:\\Program Files\\ImageMagick-7.0.8-Q16\\magick.exe" -delay 20 snap0*.png anim.gif`)