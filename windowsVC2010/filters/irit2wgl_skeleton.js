/**************************************
 *   WebGL JS skeleton for irit2wgl   *
 *   Author: Avi Kaplan.              *
 **************************************/

/*************************
 *   Enumeration Types   *
 *************************/

// Draw mode type:
var DrawMode = {
	NONE: 0x0000,
	WIREFRAME: 0x0001,
	SOLID: 0x0002,
	TEXTURE: 0x0004,

	// Special
	RENDER_TO_TEXTURE: 0x0008
};

// View angle type:
var ViewAngle = {
	ORIGINAL: 0,
	FRONT: 1,
	BACK: 2,
	RIGHT: 3,
	LEFT: 4,
	TOP: 5,
	BOTTOM: 6
};

// Projection mode type:
var ProjectionMode = {
	ORTHOGONAL: 0,
	PERSPECTIVE: 1
};

// Model object tag type:
var ModelObjectTag = {
	POLYLINE: 0,
	POLYGON: 1
};

/*****************
 *   Constants   *
 *****************/

// Dump: bkcol

// Status constants:
var STATUS_TS_MAX_ACC_LEN = 10;

var STATUS_LOG_ROW_NUM = 3;
var STATUS_LOG_COL_NUM = 0;	// Will be set after canvas width is known.

// EventMgr constants:
var KEY_B = 66;
var KEY_F = 70;
var KEY_H = 72;
var KEY_S = 83;
var KEY_Z = 90;
var KEY_LEFT = 37;
var KEY_UP = 38;
var KEY_RIGHT = 39;
var KEY_DOWN = 40;
var KEY_SHIFT = 16;
var KEY_CTRL = 17;

// Scene constants:
var SCENE_AXES_RGBA = [1, 1, 1, 0.5];

var SCENE_PICKED_INC_MATERIAL_EMISSION = [0.5, 0.5, 0.5];

// Model constants:
var MODEL_SCALE_UP_FACTOR = 0.95;
var MODEL_SCALE_DOWN_FACTOR = 1.05;

// ModelObject constants:
var MODEL_OBJECT_DEFAULT_MATERIAL_AMBIENT = [1, 1, 1];
var MODEL_OBJECT_DEFAULT_MATERIAL_DIFFUSE = [1, 1, 1];
var MODEL_OBJECT_DEFAULT_MATERIAL_SPECULAR = [1, 1, 1];
var MODEL_OBJECT_DEFAULT_MATERIAL_EMISSION = [0, 0, 0];
var MODEL_OBJECT_DEFAULT_MATERIAL_SHININESS = 8;

// Camera constants:
var VIEW_NUM = 6 + 1;

var CAMERA_VIEW_EYE = [
	[0, 0, 1],
	[0, 0, -1],
	[1, 0, 0],
	[-1, 0, 0],
	[0, 1, 0],
	[0, -1, 0]
];

var CAMERA_VIEW_AT = [
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0]
];

var CAMERA_VIEW_UP = [
	[0, 1, 0],
	[0, 1, 0],
	[0, 1, 0],
	[0, 1, 0],
	[0, 0, -1],
	[0, 0, -1]
];

var CAMERA_ORTHO_LEFT = -1;
var CAMERA_ORTHO_RIGHT = 1;
var CAMERA_ORTHO_BOTTOM = -1;
var CAMERA_ORTHO_TOP = 1;
var CAMERA_ORTHO_NEAR = -2;
var CAMERA_ORTHO_FAR = 2;

var CAMERA_PAN_STEP = 0.05;

var CAMERA_ZOOM_IN_FACTOR = 0.99;
var CAMERA_ZOOM_OUT_FACTOR = 1.01;

// Light constants:
var LIGHT_SOURCE_NUM = 8 + 1;

var LIGHT_ATTENUATION = [1, 0, 0];

var LIGHT_LAMP_RGBA = [1, 1, 1, 0.5];
var LIGHT_LAMP_MATERIAL_AMBIENT = [0, 0, 0];
var LIGHT_LAMP_MATERIAL_DIFFUSE = [0, 0, 0];
var LIGHT_LAMP_MATERIAL_SPECULAR = [1, 1, 1];
var LIGHT_LAMP_MATERIAL_EMISSION = [0.5, 0.5, 0.5];
var LIGHT_LAMP_MATERIAL_SHININESS = 64;

var LIGHT_RAY_RGBA = [1, 1, 1, 0.5];

/***************
 *   Globals   *
 ***************/

// Canvas:
var canvas = null;

// WebGL context:
var GL = null;

// Shaders:
var shaders = null;

// ControlBar:
var controlBar = null;

// Menu:
var menu = null;

// StatusLog
var statusLog = null;

// Scene:
var scene = null;

// Textures:
var textures = null;

// Picker:
var picker = null;

// Event manager:
var eventMgr = null;

/****************************
 *   Management Functions   *
 ****************************/

// Issue Warning.
function warning(msg) {
	alert("WARNING: " + msg);
}

// Report an error.
function error(msg) {
	alert("ERROR: " + msg);
	window.close();
}

// Request animation frame.
_requestAnimationFrame = (function() {
	return window.requestAnimationFrame
		|| window.webkitRequestAnimationFrame	// Chrome/Webkit
		|| window.mozRequestAnimationFrame	// Firefox
		|| window.oRequestAnimationFrame	// Opera
		|| window.msRequestAnimationFrame	// IE 10 PP2+
		|| function(callback, element) {
			window.setTimeout(callback, 1000/60);
		}
})();

// Handle keys.
function handleKeys()
{
	if (eventMgr.key.pressed[KEY_LEFT])
		scene.panCamera([CAMERA_PAN_STEP, 0, 0]);

	if (eventMgr.key.pressed[KEY_UP])
		scene.panCamera([0, -CAMERA_PAN_STEP, 0]); 

	if (eventMgr.key.pressed[KEY_RIGHT])
		scene.panCamera([-CAMERA_PAN_STEP, 0, 0]);

	if (eventMgr.key.pressed[KEY_DOWN])
		scene.panCamera([0, CAMERA_PAN_STEP, 0]);

	if (eventMgr.key.pressed[KEY_B])
		scene.panCamera([0, 0, -CAMERA_PAN_STEP]);
 
	if (eventMgr.key.pressed[KEY_F])
		scene.panCamera([0, 0, CAMERA_PAN_STEP]);

	if (eventMgr.key.pressed[KEY_Z]) {
		if (eventMgr.key.pressed[KEY_SHIFT])
			var zoomFactor = CAMERA_ZOOM_OUT_FACTOR;
		else
			var zoomFactor = CAMERA_ZOOM_IN_FACTOR;
		scene.zoomCamera(zoomFactor);
	}

	if (eventMgr.key.pressed[KEY_CTRL] && eventMgr.key.pressed[KEY_SHIFT]) {
		if (eventMgr.key.pressed[KEY_S]) {
			statusLog.showCLI = true;
		} else if (eventMgr.key.pressed[KEY_H]) {
			statusLog.showCLI = false;
		}
	}
}

// Draw.
function draw()
{
	// Clear buffers.
	GL.clear(GL.COLOR_BUFFER_BIT | GL.DEPTH_BUFFER_BIT);

	// Draw scene.
	scene.draw();

	// Finish
	GL.finish();
}

// Tick.
function tick()
{
	// TODO: invoke 'requestAnimationFrame' instead, when it will have
	//       cross-browser support.
	_requestAnimationFrame(tick);

	// Handle keys.
	handleKeys();

	// Draw.
	var preRenderTS = new Date().getTime();
	draw();
	var postRenderTS = new Date().getTime();

	// If status is hidden, no need to calculate and update.
	if (eventMgr.doc.status.view.value == "ON")
		return;

	// Calculate status.
	statusLog.TS.acc += postRenderTS - preRenderTS;
	statusLog.TS.accLen++;
	if (statusLog.TS.accLen == STATUS_TS_MAX_ACC_LEN) {
		var accSec = statusLog.TS.acc/1000;
		statusLog.FPS = STATUS_TS_MAX_ACC_LEN/(accSec);
		statusLog.TS.acc = 0;
		statusLog.TS.accLen = 0;
	}

	// Update status.
	statusLog.update();
}

// Initialize.
function init()
{
	// Get canvas.
	canvas = document.getElementById("WebGL-canvas");

	// Create WebGL context.
	try {
		GL = canvas.getContext("experimental-webgl",
			{preserveDrawingBuffer: true});
	} catch(exep) {
		error("Exception catched in getContext [" + exep.toString() + "].");
	}
	if (!GL) {
		error("Failed to create WebGL context!");
		return;
	}

	// Initialize the shaders.
	shaders = new Shaders();
	shaders.init();

	// Initialize the control bar.
	controlBar = new ControlBar();
	controlBar.init();

	// Initialize the menu.
	menu = new Menu();
	menu.init();

	// Initialize the status log.
	statusLog = new StatusLog();
	statusLog.init();

	// Initialize the textures.
	textures = new Textures();
	textures.init();

	// Initialize the scene.
	scene = new Scene();
	scene.init();

	// Initialize the picker.
	picker = new Picker();
	picker.init();

	// Initialize the event manager.
	eventMgr = new EventMgr();
	eventMgr.init();

	// Set viewport.
	GL.viewport(0, 0, canvas.width, canvas.height);

	// Clear background color.
	GL.clearColor(BK_RGB[0]/255, BK_RGB[1]/255, BK_RGB[2]/255, 1);

	// Clear depth.
	GL.clearDepth(1);

	// Set polygon offset.
	GL.polygonOffset(1, 0);
	GL.enable(GL.POLYGON_OFFSET_FILL);

	// Run.
	tick();
}

// Destroy.
function destroy()
{
	// Destroy the event manager.
	eventMgr.destroy();

	// Destroy the picker.
	picker.destroy();

	// Destroy the scene.
	scene.destroy();

	// Destroy the textures.
	textures.destroy();

	// Destroy the status log.
	statusLog.destroy();

	// Destroy the menu.
	menu.destroy();

	// Destroy the control bar.
	controlBar.destroy();

	// Destroy the shaders.
	shaders.destroy();
}

/*********************
 *   Class Shaders   *
 *********************/

// Constructor.
Shaders = function() {
	this.program = null;
	this.loc = new Object();
};

// Initialize the shaders.
Shaders.prototype.init = function() {
	// Create and compile the vertex shader.
	var vertexShader = GL.createShader(GL.VERTEX_SHADER);
	GL.shaderSource(vertexShader, this.getShaderSource("vshader"));
	GL.compileShader(vertexShader);
	if (!GL.getShaderParameter(vertexShader, GL.COMPILE_STATUS)) {
		error("Error during vertex shader compilation:\n"
			+ GL.getShaderInfoLog(vertexShader));
		return;
	}

	// Create and compile the fragment shader.
	var fragmentShader = GL.createShader(GL.FRAGMENT_SHADER);
	GL.shaderSource(fragmentShader, this.getShaderSource("fshader"));
	GL.compileShader(fragmentShader);
	if (!GL.getShaderParameter(fragmentShader, GL.COMPILE_STATUS)) {
		error("Error during fragment shader compilation:\n"
			+ GL.getShaderInfoLog(fragmentShader));
		return;
	}

	// Create the shaders program and link it.
	this.program = GL.createProgram();
	GL.attachShader(this.program, fragmentShader);
	GL.attachShader(this.program, vertexShader);
	GL.linkProgram(this.program);
	if (!GL.getProgramParameter(this.program, GL.LINK_STATUS)) {
		error("Error during program linking:\n"
			+ GL.getProgramInfoLog(this.program));
		return;
	}

	// Validate the shaders program.
	GL.validateProgram(this.program);
	if (!GL.getProgramParameter(this.program, GL.VALIDATE_STATUS)) {
		error("Error during program validation:\n"
			+ GL.getProgramInfoLog(this.program));
		return;
	}

	// Use the shaders program.
	GL.useProgram(this.program);

	// Set locations.
	this.setLocations();
};

// Destroy the shaders.
Shaders.prototype.destroy = function() {
	// Nothing to do here ...
};

// Get shader source code.
Shaders.prototype.getShaderSource = function(shaderId) {
	var shaderScript = document.getElementById(shaderId);
	if (!shaderScript)
		return null;

	var sourceStr = "";
	var node = shaderScript.firstChild;
	while (node) {
		if (node.nodeType == Node.TEXT_NODE) 
			sourceStr += node.textContent;
		node = node.nextSibling;
	}

	return sourceStr;
};

// Get attribute location.
Shaders.prototype.getAttribLocation = function(attr) {
	var loc = GL.getAttribLocation(this.program, attr);
	if (loc == -1)
		error("Error during " + attr + " address retrival.");
	return loc;
};

// Get uniform location.
Shaders.prototype.getUniformLocation = function(attr) {
	var loc = GL.getUniformLocation(this.program, attr);
	if (loc == -1)
		error("Error during " + attr + " address retrival.");
	return loc;
};

// Set locations.
Shaders.prototype.setLocations = function() {
	// Vertex position.
	this.loc.aVertexPosition = this.getAttribLocation('aVertexPosition'); 

	// Vertex normal.
	this.loc.aVertexNormal = this.getAttribLocation('aVertexNormal');

	// Vertex texture UV.
	this.loc.aVertexTextureUV = this.getAttribLocation('aVertexTextureUV');

	// Model matrix.
	this.loc.uMMat = this.getUniformLocation('uMMat');

	// View matrix.
	this.loc.uVMat = this.getUniformLocation('uVMat');

	// Projection matrix.
	this.loc.uPMat = this.getUniformLocation('uPMat');

	// Normal matrix.
	this.loc.uNMat = this.getUniformLocation('uNMat');

	// COP (center of projeciton).
	this.loc.uCOP = this.getUniformLocation('uCOP');

	// RGBA.
	this.loc.uRGBA = this.getUniformLocation('uRGBA');

	// Apply texture?
	this.loc.uApplyTexture = this.getUniformLocation('uApplyTexture');

	// Sampler.
	this.loc.uSampler = this.getUniformLocation('uSampler');

	// Material.
	this.loc.uMaterial = new Object();
	this.loc.uMaterial.ambient = this.getUniformLocation('uMaterial.ambient');
	this.loc.uMaterial.diffuse = this.getUniformLocation('uMaterial.diffuse');
	this.loc.uMaterial.specular
		= this.getUniformLocation('uMaterial.specular');
	this.loc.uMaterial.emission
		= this.getUniformLocation('uMaterial.emission');
	this.loc.uMaterial.shininess
		= this.getUniformLocation('uMaterial.shininess');

	// Use light?
	this.loc.uUseLight = this.getUniformLocation('uUseLight');

	// Lights.
	this.loc.uLights = new Array();
	for (var i = 0; i < LIGHT_SOURCE_NUM; i++) {
		this.loc.uLights[i] = new Object();
		this.loc.uLights[i].type
			= this.getUniformLocation('uLights[' + i + '].type');
		this.loc.uLights[i].TL
			= this.getUniformLocation('uLights[' + i + '].TL');
		this.loc.uLights[i].pos
			= this.getUniformLocation('uLights[' + i + '].pos');
		this.loc.uLights[i].dir
			= this.getUniformLocation('uLights[' + i + '].dir');
		this.loc.uLights[i].attK
			= this.getUniformLocation('uLights[' + i + '].attK');
		this.loc.uLights[i].ambient
			= this.getUniformLocation('uLights[' + i + '].ambient');
		this.loc.uLights[i].diffuse
			= this.getUniformLocation('uLights[' + i + '].diffuse');
		this.loc.uLights[i].specular
			= this.getUniformLocation('uLights[' + i + '].specular');
	}
};

/***********************
 *   Class ControlBar   *
 ***********************/

// Constructor.
ControlBar = function() {
};

// Initialize the control bar.
ControlBar.prototype.init = function() {
	// Set parameters.
	this.setParams();

	// Show/Hide control bar?
	if (this.show)
		document.getElementById("control-data").style.display = "block";
	else
		document.getElementById("control-data").style.display = "none";
};

// Destroy the control bar.
ControlBar.prototype.destroy = function() {
	// Nothing to do here ...
};

// Dump: ctrl

/******************
 *   Class Menu   *
 ******************/

 // Constructor.
Menu = function() {
	this.scene = new Object();
	this.model = new Object();
	this.camera = new Object();
	this.light = new Object();
};

// Initialize the menu.
Menu.prototype.init = function() {
	// Set parameters.
	this.setParams();
};

// Destroy the menu.
Menu.prototype.destroy = function() {
	// Nothing to do here ...
};

// Dump: menu

// Update menu parameters.
Menu.prototype.updateParams = function() {
	// Scene menu parameters.
	this.scene.showAxes = scene.worldAxes.visible;
	this.scene.enableDepthTest = scene.enableDepthTest;
	this.scene.enablePicking = scene.enablePicking;

	// Model menu parameters.
	this.model.showAxes = scene.modelAxes.visible;
	this.model.drawMode = scene.drawMode;
	this.model.worldTrans = scene.model.worldTrans;

	// Camera menu parameters.
	this.camera.viewAngle = scene.camera.viewAngle;
	this.camera.projMode = scene.camera.projMode;

	// Light menu parameters.
	this.light.ambient = scene.lights[0].ambient;
	this.light.ambientI = scene.lights[0].ambientI;
};

Menu.prototype.getParamsInUsageFormat = function() {
	this.updateParams();

	var usageFormat =
		(this.scene.showAxes ? ' -W' : '')
		+ (this.scene.enableDepthTest ? '' : ' -D')
		+ (this.scene.enablePicking ? ' -P' : '')
		+ (this.model.showAxes ? ' -M' : '')
		+ ' -d ' + this.model.drawMode
		+ (this.model.worldTrans ? '' : ' -T')
		+ ' -v ' + this.camera.viewAngle
		+ ' -p ' + this.camera.projMode
		+ ' -a ' + this.light.ambient[0]
		+ ' ' + this.light.ambient[1]
		+ ' ' + this.light.ambient[2];

		return usageFormat;
};

/***********************
 *   Class StatusLog   *
 ***********************/

// Constructor.
StatusLog = function() {
	// Number of polygons.
	this.NOP = 0;

	// Frames per second.
	this.FPS = 0;

	// Resolution of tesselation.
	this.resolution = 0;

	// Show CLI menu parameters?
	this.showCLI = false;

	this.TS = new Object();

	this.doc = new Object();
};

// Initialize the status log.
StatusLog.prototype.init = function() {
	this.doc.log = document.getElementById("status-log");
	this.doc.log.rows = STATUS_LOG_ROW_NUM;
	STATUS_LOG_COL_NUM = Math.floor(canvas.width/8.2);
	this.doc.log.cols = STATUS_LOG_COL_NUM;

	this.TS.acc = 0;
	this.TS.accLen = 0;

	// Set parameters.
	this.setParams();
};

// Destroy the status log.
StatusLog.prototype.destroy = function() {
	// Nothing to do here ...
};

// Dump: log

// Update the status log.
StatusLog.prototype.update = function() {
	var row1 = 'NOP: ' + this.NOP.toString()
		+ '   Resolution: ' + this.resolution.toString()
		+ '   FPS: ' + this.FPS.toFixed(5)
		+ '\n'

	var row2 = picker.objName + '\n';

	var row3 = '';

	if (this.showCLI) {
		row3 += 'CLI: \"'
			+ (controlBar.show ? '': ' -C ')
			+ ' -w ' + canvas.width + ' -h ' + canvas.height
			+ ' -b ' + BK_RGB[0] + ' ' + BK_RGB[1] + ' ' + BK_RGB[2]
			+ menu.getParamsInUsageFormat() + '\"';
	}

	this.doc.log.cols = Math.max(row1.length, row2.length, row3.length,
		STATUS_LOG_COL_NUM);

	this.doc.log.value = row1 + row2 + row3;
};

/**********************
 *   Class EventMgr   *
 **********************/
 
// Constructor.
EventMgr = function() {
	this.mouse = new Object();
	this.key = new Object();

	this.doc = new Object();
};

// Initialize the event manager.
EventMgr.prototype.init = function() {
	this.mouse.rbDown = false;
	this.mouse.mbDown = false;
	this.mouse.lbDown = false;

	this.key.pressed = {};

	// Set the event handlers.
	this.setEventHandlers();

	// Set light source menu.
	this.setLightSourceMenu();
};

// Destroy the event manager.
EventMgr.prototype.destroy = function() {
	// Nothing to do here ...
};

// Set the event handlers.
EventMgr.prototype.setEventHandlers = function() {
	this.setControlEventHandlers();
	this.setMouseEventHandlers();
	this.setKeyboardEventHandlers();
	this.setSceneEventHandlers();
	this.setModelEventHandlers();
	this.setCameraEventHandlers();
	this.setLightEventHandlers();
};

// Set control event handlers.
EventMgr.prototype.setControlEventHandlers = function() {
	this.doc.menu = new Object();
	this.doc.menu.view = document.getElementById("menu-view");
	this.doc.menu.view.onclick = handleMenuViewOnClick;
	this.doc.menu.view.value = "ON";
	document.getElementById("menu-data").style.display = "none";

	this.doc.status = new Object();
	this.doc.status.view = document.getElementById("status-view");
	this.doc.status.view.onclick = handleStatusViewOnClick;
	this.doc.status.view.value = "ON";
	document.getElementById("status-data").style.display = "none";
};

// Set mouse event handlers.
EventMgr.prototype.setMouseEventHandlers = function() {
	document.onmouseup = handleMouseUp;
	canvas.onmousemove = handleMouseMove;
	canvas.onmousedown = handleMouseDown;
	canvas.onmousewheel = handleMouseWheel;
};

// Set keyboard event handlers.
EventMgr.prototype.setKeyboardEventHandlers = function() {
	document.onkeydown = handleKeyDown;
	document.onkeyup = handleKeyUp;
};

// Set scene event handlers.
EventMgr.prototype.setSceneEventHandlers = function() {
	this.doc.scene = new Object();

	// Axes.
	this.doc.scene.axes = document.getElementById("scene-axes");
	this.doc.scene.axes.onclick = handleSceneAxesOnClick;
	this.doc.scene.axes.checked = menu.scene.showAxes;

	// Depth.
	this.doc.scene.depth = document.getElementById("scene-depth");
	this.doc.scene.depth.onclick = handleSceneDepthOnClick;
	this.doc.scene.depth.checked = menu.scene.enableDepthTest;

	// Picking.
	this.doc.scene.picking = document.getElementById("scene-picking");
	this.doc.scene.picking.onclick = handleScenePickingOnClick;
	this.doc.scene.picking.checked = menu.scene.enablePicking;
};

// Set model event handlers.
EventMgr.prototype.setModelEventHandlers = function() {
	this.doc.model = new Object();

	// Axes.
	this.doc.model.axes = document.getElementById("model-axes");
	this.doc.model.axes.onclick = handleModelAxesOnClick;
	this.doc.model.axes.checked = menu.model.showAxes;

	// Draw mode.
	this.doc.model.drawmode = new Object();

	this.doc.model.drawmode.wireframe
		= document.getElementById("model-dm-wire");
	this.doc.model.drawmode.wireframe.onclick = handleModelDrawModeOnClick;
	this.doc.model.drawmode.wireframe.checked
		= (menu.model.drawMode & DrawMode.WIREFRAME);

	this.doc.model.drawmode.solid
		= document.getElementById("model-dm-solid");
	this.doc.model.drawmode.solid.onclick = handleModelDrawModeOnClick;
	this.doc.model.drawmode.solid.checked
		= (menu.model.drawMode & DrawMode.SOLID);

	this.doc.model.drawmode.texture
		= document.getElementById("model-dm-texture");
	this.doc.model.drawmode.texture.onclick = handleModelDrawModeOnClick;
	this.doc.model.drawmode.texture.checked
		= (menu.model.drawMode & DrawMode.TEXTURE);

	// Transformation.
	this.doc.model.trans = new Object();

	this.doc.model.trans.world = document.getElementById("model-trans-world");
	this.doc.model.trans.world.onclick = handleModelTransformationOnClick;
	this.doc.model.trans.world.checked = menu.model.worldTrans;

	this.doc.model.trans.object = document.getElementById("model-trans-object");
	this.doc.model.trans.object.onclick = handleModelTransformationOnClick;
	this.doc.model.trans.object.checked = !(menu.model.worldTrans);
};

// Set camera event handlers.
EventMgr.prototype.setCameraEventHandlers = function() {
	this.doc.camera = new Object();

	// View.
	this.doc.camera.view = new Object();

	this.doc.camera.view.which = document.getElementById("camera-view-which");
	this.doc.camera.view.which.onchange = handleCameraViewWhichOnChange;
	this.doc.camera.view.which.value = menu.camera.viewAngle.toString();

	// Projection.
	this.doc.camera.projmode = new Object();

	this.doc.camera.projmode.orthogonal
		= document.getElementById("camera-pm-orthogonal");
	this.doc.camera.projmode.orthogonal.onclick = handleCameraProjModeOnClick;
	this.doc.camera.projmode.orthogonal.checked
		= (menu.camera.projMode == ProjectionMode.ORTHOGONAL);

	this.doc.camera.projmode.perspective
		= document.getElementById("camera-pm-perspective");
	this.doc.camera.projmode.perspective.onclick = handleCameraProjModeOnClick;
	this.doc.camera.projmode.perspective.checked
		= (menu.camera.projMode == ProjectionMode.PERSPECTIVE);
};

// Set light event handlers.
EventMgr.prototype.setLightEventHandlers = function() {
	this.doc.light = new Object();

	// Light ambient.
	this.doc.light.ambient = document.getElementById("light-ambient");
	this.doc.light.ambient.onchange = handleLightAmbientOnChange;
	this.doc.light.ambient.value = scene.lights[0].ambientI*255;

	// Light source.
	this.doc.light.src = new Object();
	
	this.doc.light.src.which = document.getElementById("light-src-which");
	this.doc.light.src.which.onchange = handleLightSrcWhichOnChange;
	this.doc.light.src.which.value = "1";

	this.doc.light.src.visible = document.getElementById("light-src-visible");
	this.doc.light.src.visible.onclick = handleLightSrcVisibleOnClick;
	this.doc.light.src.visible.checked = false;

	// Light type.
	this.doc.light.type = new Object();

	this.doc.light.type.point = document.getElementById("light-type-point");
	this.doc.light.type.point.onclick = handleLightTypeOnClick;
	this.doc.light.type.point.checked
		= (scene.lights[1].type == LightType.POINT);

	this.doc.light.type.directional
		= document.getElementById("light-type-directional");
	this.doc.light.type.directional.onclick = handleLightTypeOnClick;
	this.doc.light.type.directional.checked = !(this.doc.light.type.point.checked);
	
	// Light TL (2 lights).
	this.doc.light.TL = document.getElementById("light-TL");
	this.doc.light.TL.onclick = handleLightTLOnClick;
	this.doc.light.TL.checked = scene.lights[1].TL;

	// Light position.
	this.doc.light.pos = new Object();

	this.doc.light.pos.set = document.getElementById("light-pos-set");
	this.doc.light.pos.set.onclick = handleLightPosSetOnClick;

	eventMgr.doc.light.pos.x = document.getElementById("light-pos-X");
	eventMgr.doc.light.pos.y = document.getElementById("light-pos-Y");
	eventMgr.doc.light.pos.z = document.getElementById("light-pos-Z");

	// Light direciton.
	this.doc.light.dir = new Object();

	this.doc.light.dir.set = document.getElementById("light-dir-set");
	this.doc.light.dir.set.onclick = handleLightDirSetOnClick;

	eventMgr.doc.light.dir.x = document.getElementById("light-dir-X");
	eventMgr.doc.light.dir.y = document.getElementById("light-dir-Y");
	eventMgr.doc.light.dir.z = document.getElementById("light-dir-Z");

	// Light intensity.
	this.doc.light.intens = new Object();

	this.doc.light.intens.diffuse
		= document.getElementById("light-intens-diffuse");
	this.doc.light.intens.diffuse.onchange = handleLightIntensDiffuseOnChange;
	this.doc.light.intens.specular
		= document.getElementById("light-intens-specular");
	this.doc.light.intens.specular.onchange = handleLightIntensSpecularOnChange;
};

// Set light source menu.
EventMgr.prototype.setLightSourceMenu = function() {
	var lightId = parseInt(eventMgr.doc.light.src.which.value);

	this.doc.light.intens.diffuse.value = scene.lights[lightId].diffuseI*255;
	this.doc.light.intens.specular.value = scene.lights[lightId].specularI*255;

	switch (scene.lights[lightId].type) {
	case LightType.POINT:
		document.getElementById("light-pos").style.display = "block";
		document.getElementById("light-dir").style.display = "none";

		eventMgr.doc.light.pos.x.value = scene.lights[lightId].pos[0];
		eventMgr.doc.light.pos.y.value = scene.lights[lightId].pos[1];
		eventMgr.doc.light.pos.z.value = scene.lights[lightId].pos[2];
		break;

	case LightType.DIRECTIONAL:
		document.getElementById("light-pos").style.display = "none";
		document.getElementById("light-dir").style.display = "block";

		eventMgr.doc.light.dir.x.value = scene.lights[lightId].dir[0];
		eventMgr.doc.light.dir.y.value = scene.lights[lightId].dir[1];
		eventMgr.doc.light.dir.z.value = scene.lights[lightId].dir[2];
		break;
	}
};

// Handle menu view click events.
handleMenuViewOnClick = function(event)
{
	if (eventMgr.doc.menu.view.value == "ON") {
		document.getElementById("menu-data").style.display = "block";
		eventMgr.doc.menu.view.value = "OFF";
	} else {
		document.getElementById("menu-data").style.display = "none"; 
		eventMgr.doc.menu.view.value = "ON";
	}
}

// Handle status view click events.
handleStatusViewOnClick = function(event)
{
	if (eventMgr.doc.status.view.value == "ON") {
		document.getElementById("status-data").style.display = "block";
		eventMgr.doc.status.view.value = "OFF";
	} else {
		document.getElementById("status-data").style.display = "none"; 
		eventMgr.doc.status.view.value = "ON";
	}
}

// Handle mouse down events.
handleMouseDown = function(event)
{
	var currX = event.offsetX ? event.offsetX : event.layerX;
	var currY = event.offsetY ? event.offsetY : event.layerY;

	if (event.button == 0) {
		eventMgr.mouse.lbDown = true;
	} if (event.button == 1)
		eventMgr.mouse.mbDown = true;
	else if (event.button == 2) {
		eventMgr.mouse.rbDown = true;
		if (scene.enablePicking)
			picker.pick(currX, canvas.height - currY);
	}

	eventMgr.mouse.lastX = currX;
	eventMgr.mouse.lastY = currY;

	// Prevent default mouse event handling.
	if (event.preventDefault)
		event.preventDefault();
	else
		event.returnValue = false;
    return false;
}

// Handle mouse up events.
function handleMouseUp(event)
{
	if (event.button == 0)
		eventMgr.mouse.lbDown = false;
	if (event.button == 1)
		eventMgr.mouse.mbDown = false;
	else if (event.button == 2)
		eventMgr.mouse.rbDown = false;
}

// Handle mouse move events.
function handleMouseMove(event)
{
	var currX = event.offsetX ? event.offsetX : event.layerX;
	var currY = event.offsetY ? event.offsetY : event.layerY;

	var deltaX = currX - eventMgr.mouse.lastX;
	var deltaY = currY - eventMgr.mouse.lastY;

	var delta = [2*deltaX/canvas.width, -2*deltaY/canvas.height, 0];

	if (eventMgr.mouse.lbDown) {
		if (eventMgr.key.pressed[17])	// ctrl
			scene.spinCamera([-delta[1], delta[0], 0]);
		else
			scene.rotateModel(delta);
	}

	if (eventMgr.mouse.mbDown)
		scene.scaleModel(delta);

	if (eventMgr.mouse.rbDown) {
		if (eventMgr.key.pressed[17])	// ctrl
			scene.panCamera(delta);
		else
			scene.translateModel(delta);
	}

	eventMgr.mouse.lastX = currX;
	eventMgr.mouse.lastY = currY;
}

// Handle mouse wheel events.
function handleMouseWheel(event)
{
	// TODO
}

// Handle key down events.
function handleKeyDown(event)
{
	var keyVal = (event.which) ? event.which : event.keyCode;
	eventMgr.key.pressed[keyVal] = true;
}

// Handle key up events.
function handleKeyUp(event)
{
	var keyVal = (event.which) ? event.which : event.keyCode;
	eventMgr.key.pressed[keyVal] = false;
}

// Handle scene axes click events.
function handleSceneAxesOnClick(event)
{
	scene.worldAxes.visible = eventMgr.doc.scene.axes.checked;
}

// Handle scene depth click events.
function handleSceneDepthOnClick(event)
{
	if (eventMgr.doc.scene.depth.checked)
		scene.enableDepthTest = true;
	else
		scene.enableDepthTest = false;
}

// Handle scene picking click events.
function handleScenePickingOnClick(event)
{
	if (eventMgr.doc.scene.picking.checked)
		scene.enablePicking = true;
	else {
		scene.enablePicking = false;
		picker.applyPick(-1);
	}
}

// Handle model axes click events.
function handleModelAxesOnClick(event)
{
	scene.modelAxes.visible = eventMgr.doc.model.axes.checked;
}

// Handle model draw-mode click events.
function handleModelDrawModeOnClick(event)
{
	if (eventMgr.doc.model.drawmode.wireframe.checked)
		scene.drawMode |= DrawMode.WIREFRAME;
	else
		scene.drawMode &= ~DrawMode.WIREFRAME;

	if (eventMgr.doc.model.drawmode.solid.checked)
		scene.drawMode |= DrawMode.SOLID;
	else
		scene.drawMode &= ~DrawMode.SOLID;

	if (eventMgr.doc.model.drawmode.texture.checked)
		scene.drawMode |= DrawMode.TEXTURE;
	else
		scene.drawMode &= ~DrawMode.TEXTURE;
}

// Handle model transformation click events.
function handleModelTransformationOnClick(event)
{
	if (eventMgr.doc.model.trans.world.checked)
		scene.model.worldTrans = true;
	else if (eventMgr.doc.model.trans.object.checked)
		scene.model.worldTrans = false;
}

// Handle camera view change events.
function handleCameraViewWhichOnChange(event)
{
	scene.camera.viewAngle = parseInt(eventMgr.doc.camera.view.which.value);
	scene.camera.setView();
}

// Handle camera projection click events.
function handleCameraProjModeOnClick(event)
{
	if (eventMgr.doc.camera.projmode.orthogonal.checked)
		scene.camera.projMode = ProjectionMode.ORTHOGONAL;
	else if (eventMgr.doc.camera.projmode.perspective.checked)
		scene.camera.projMode = ProjectionMode.PERSPECTIVE;

	scene.camera.setProjection();
}

// Handle light ambient change events.
function handleLightAmbientOnChange(event)
{
	scene.lights[0].ambientI
		= parseInt(eventMgr.doc.light.ambient.value)/255;
}

// Handle light source which change events.
function handleLightSrcWhichOnChange(event)
{
	var lightId = parseInt(eventMgr.doc.light.src.which.value);

	eventMgr.doc.light.src.visible.checked = scene.lights[lightId].lamp.visible;

	switch (scene.lights[lightId].type) {
	case LightType.POINT:
		eventMgr.doc.light.type.point.checked = true;
		eventMgr.doc.light.type.directional.checked = false;
		break;

	case LightType.DIRECTIONAL:
	eventMgr.doc.light.type.directional.checked = false;
		eventMgr.doc.light.type.directional.checked = true;
		break;
	}

	eventMgr.doc.light.TL.checked = scene.lights[lightId].TL;

	eventMgr.setLightSourceMenu();
}

// Handle light source visible click events.
function handleLightSrcVisibleOnClick(event)
{
	var lightId = parseInt(eventMgr.doc.light.src.which.value);

	scene.lights[lightId].lamp.visible = eventMgr.doc.light.src.visible.checked;
	scene.lights[lightId].ray.visible = eventMgr.doc.light.src.visible.checked;
}

// Handle light type click events.
function handleLightTypeOnClick(event)
{
	var lightId = parseInt(eventMgr.doc.light.src.which.value);

	if (eventMgr.doc.light.type.point.checked) {
		scene.lights[lightId].type = LightType.POINT;
	} else if (eventMgr.doc.light.type.directional.checked) {
		scene.lights[lightId].type = LightType.DIRECTIONAL;
	}

	eventMgr.setLightSourceMenu();
}

// Handle light TL click events.
function handleLightTLOnClick(event)
{
	var lightId = parseInt(eventMgr.doc.light.src.which.value);

	scene.lights[lightId].TL = eventMgr.doc.light.TL.checked;
}

// Handle light position set click events.
function handleLightPosSetOnClick(event)
{
	posX = parseFloat(eventMgr.doc.light.pos.x.value);
	posY = parseFloat(eventMgr.doc.light.pos.y.value);
	posZ = parseFloat(eventMgr.doc.light.pos.z.value);
	if (isNaN(posX) | isNaN(posY) | isNaN(posZ)) {
		alert("Invalid \'position\' value!");
		return;
	}

	var lightId = parseInt(eventMgr.doc.light.src.which.value);
	scene.lights[lightId].setPos([posX, posY, posZ]);
}

// Handle light direction set click events.
function handleLightDirSetOnClick(event)
{
	dirX = parseFloat(eventMgr.doc.light.dir.x.value);
	dirY = parseFloat(eventMgr.doc.light.dir.y.value);
	dirZ = parseFloat(eventMgr.doc.light.dir.z.value);
	if (isNaN(dirX) | isNaN(dirY) | isNaN(dirZ)) {
		alert("Invalid \'direction\' value!");
		return;
	}

	var lightId = parseInt(eventMgr.doc.light.src.which.value);
	scene.lights[lightId].setDir([dirX, dirY, dirZ]);
}

// Handle light intensity diffuse change events.
function handleLightIntensDiffuseOnChange(event)
{
	var lightId = parseInt(eventMgr.doc.light.src.which.value);
	scene.lights[lightId].diffuseI
		= parseInt(eventMgr.doc.light.intens.diffuse.value)/255;
}

// Handle light intensity specular change events.
function handleLightIntensSpecularOnChange(event)
{
	var lightId = parseInt(eventMgr.doc.light.src.which.value);
	scene.lights[lightId].specularI
		= parseInt(eventMgr.doc.light.intens.specular.value)/255;
}

/*******************
*   Class Picker   *
 *******************/

// Constructor.
Picker = function() {
	this.frameBuffer = null;
	this.renderBuffer = null;
	this.texture = null;

	this.objId = -1;
	this.objName = '';
};

// Initialize the picker.
Picker.prototype.init = function() {
	this.frameBuffer = GL.createFramebuffer();
	GL.bindFramebuffer(GL.FRAMEBUFFER, this.frameBuffer);

	this.texture = GL.createTexture();
	GL.bindTexture(GL.TEXTURE_2D, this.texture);
	GL.texParameteri(GL.TEXTURE_2D, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
	GL.texParameteri(GL.TEXTURE_2D, GL.TEXTURE_MIN_FILTER,
		GL.LINEAR_MIPMAP_NEAREST);
	GL.generateMipmap(GL.TEXTURE_2D);

	GL.texImage2D(GL.TEXTURE_2D, 0, GL.RGBA, canvas.width, canvas.height, 0,
		GL.RGBA, GL.UNSIGNED_BYTE, null);

	this.renderBuffer = GL.createRenderbuffer();
		GL.bindRenderbuffer(GL.RENDERBUFFER, this.renderBuffer);
	GL.renderbufferStorage(GL.RENDERBUFFER, GL.DEPTH_COMPONENT16, canvas.width,
		canvas.height);

	GL.framebufferTexture2D(GL.FRAMEBUFFER, GL.COLOR_ATTACHMENT0, GL.TEXTURE_2D,
		this.texture, 0); 
	GL.framebufferRenderbuffer(GL.FRAMEBUFFER, GL.DEPTH_ATTACHMENT,
		GL.RENDERBUFFER, this.renderBuffer);

	GL.bindTexture(GL.TEXTURE_2D, null); 
	GL.bindRenderbuffer(GL.RENDERBUFFER, null);
	GL.bindFramebuffer(GL.FRAMEBUFFER, null);
};

// Destroy the picker.
Picker.prototype.destroy = function() {
	deleteRenderbuffer(this.renderBuffer);
	deleteTexture(this.texture);
	deleteFramebuffer(this.frameBuffer);
};

// Pick.
Picker.prototype.pick = function(x, y) {
	GL.bindFramebuffer(GL.FRAMEBUFFER, this.frameBuffer);

	scene.drawMode |= DrawMode.RENDER_TO_TEXTURE;
	draw();
	scene.drawMode &= ~DrawMode.RENDER_TO_TEXTURE;

	var RGBA = new Uint8Array(4);
	GL.readPixels(x, y, 1, 1, GL.RGBA, GL.UNSIGNED_BYTE, RGBA);

	var objId = picker.RGBA2ID(RGBA);
	this.applyPick(objId);

	GL.bindFramebuffer(GL.FRAMEBUFFER, null);
};

// Apply pick.
Picker.prototype.applyPick = function(objId) {
	if (this.objId >= 0) {
		scene.model.objs[this.objId].material.emission = Vec3.sub(
			scene.model.objs[this.objId].material.emission,
			SCENE_PICKED_INC_MATERIAL_EMISSION
		);
	}

	if (objId >= 0) {
		scene.model.objs[objId].material.emission = Vec3.add(
			scene.model.objs[objId].material.emission,
			SCENE_PICKED_INC_MATERIAL_EMISSION
		);
		this.objName = scene.model.objs[objId].name;
	} else {
		this.objName = '';
	}

	this.objId = objId;
};

// ID -> RGBA.
Picker.prototype.ID2RGBA = function(id) {
	var RGBA = [0, 0, 0, 0];

	RGBA[0] = Math.floor( (id + 1) / (256*256) );
	RGBA[1] = Math.floor( ((id + 1) % (256*256) ) / 256);
	RGBA[2] = ((id + 1) % (256*256)) % 256;
	RGBA[3] = 255;

	return RGBA;
};

// RGBA -> ID.
Picker.prototype.RGBA2ID = function(RGBA) {
	id = RGBA[0]*256*256 + RGBA[1]*256 + RGBA[2] - 1;
	return id;
};

/**********************
 *   Class Textures   *
 **********************/

// Constructor.
Textures = function() {
	this.data = new Array();
};

// Initialize the textures.
Textures.prototype.init = function() {
	this.setData();
};

// Destroy the textures.
Textures.prototype.destroy = function() {
	for (var i = 0; i < this.data.length; i++)
		GL.deleteTexture(this.data[i]);
};

Textures.prototype.add = function(imageSrc) {
	var ext = imageSrc.split('.')[imageSrc.split('.').length - 1];

	var texture = GL.createTexture();
	texture.image = new Image();
	texture.image.onload = function() {
		handleLoadedTexture(texture)
	};
	texture.image.src = imageSrc;
	this.data.push(texture);
};

// Handle loaded texture.
function handleLoadedTexture(texture) {
	GL.bindTexture(GL.TEXTURE_2D, texture);
	GL.pixelStorei(GL.UNPACK_FLIP_Y_WEBGL, true);
	GL.texImage2D(GL.TEXTURE_2D, 0, GL.RGBA, GL.RGBA, GL.UNSIGNED_BYTE,
		texture.image);
	GL.texParameteri(GL.TEXTURE_2D, GL.TEXTURE_MAG_FILTER, GL.NEAREST);
	GL.texParameteri(GL.TEXTURE_2D, GL.TEXTURE_MIN_FILTER,
		GL.LINEAR_MIPMAP_NEAREST);
	GL.generateMipmap(GL.TEXTURE_2D);
	GL.bindTexture(GL.TEXTURE_2D, null);
}

// Dump: texture

/*******************
 *   Class Scene   *
 *******************/

// Constructor.
Scene = function() {
	this.enableDepthTest = false;
	this.enablePicking = false;

	this.model = null;

	this.drawMode = DrawMode.NONE;

	this.camera = null;

	this.lights = new Array();

	this.worldAxes = null;
	this.modelAxes = null;
};

// Initialize the scene.
Scene.prototype.init = function() {
	// Scene:
	this.enableDepthTest = menu.scene.enableDepthTest;
	this.enablePicking = menu.scene.enablePicking;

	// Model:
	this.setModel();

	// Camera:
	this.setCamera();

	// Light:
	this.setLight();

	// World axes:
	this.worldAxes = new Axes();
	this.worldAxes.visible = menu.scene.showAxes;

	// Model axes:
	this.modelAxes = new Axes();
	this.modelAxes.visible = menu.model.showAxes;
};

// Destroy the scene.
Scene.prototype.destroy = function() {
	// Model axes:
	this.modelAxes.destroyBuffers();

	// World axes:
	this.worldAxes.destroyBuffers();

	// Light:
	for (var i = 1; i < LIGHT_SOURCE_NUM; i++) {
		this.lights[i].lamp.destroyBuffers();
		if (this.lights[i].ray != null)
			this.lights[i].ray.destroyBuffers();
	}

	// Camera:
	// Nothing to do here ...

	// Model:
	this.model.destroyBuffers();
};

// Set the model.
Scene.prototype.setModel = function() {
	this.model = new Model();

	this.model.setData();
	this.model.initBuffers();

	this.drawMode = menu.model.drawMode;
	this.model.worldTrans = menu.model.worldTrans;
};

// Set the camera.
Scene.prototype.setCamera = function() {
	this.camera = new Camera();
	this.camera.initMatrices();
	this.camera.setView();
	this.camera.setProjection();
};

// Set the light.
Scene.prototype.setLight = function() {
	this.setGlobalAmbientLight();
	this.setLightSources();
};

// Set global ambient light.
Scene.prototype.setGlobalAmbientLight = function() {
	this.lights[0] = new Light(0);
	this.lights[0].type = LightType.AMBIENT;
	this.lights[0].ambient = menu.light.ambient;
	this.lights[0].ambientI = Math.max(
		this.lights[0].ambient[0],
		this.lights[0].ambient[1],
		this.lights[0].ambient[2]
	);
};

// Dump: light

// Draw scene.
Scene.prototype.draw = function() {
	// Set camera.
	this.camera.set();

	// Is render to texture mode on?
	if (this.drawMode & DrawMode.RENDER_TO_TEXTURE) {
		GL.enable(GL.DEPTH_TEST);
		GL.disable(GL.BLEND);
		this.model.draw(DrawMode.RENDER_TO_TEXTURE);
		return;
	}

	// Set lights.
	for (var i = 0; i < LIGHT_SOURCE_NUM; i++)
		this.lights[i].set();

	// Enable depth test?
	if (scene.enableDepthTest)
		GL.enable(GL.DEPTH_TEST);
	else
		GL.disable(GL.DEPTH_TEST);

	// Enable blending.
	GL.enable(GL.BLEND);
	GL.blendFunc(GL.SRC_ALPHA, GL.ONE_MINUS_SRC_ALPHA);

	// Draw texture model?
	if (this.drawMode & DrawMode.TEXTURE)
		this.model.draw(DrawMode.TEXTURE);

	// Draw solid model?
	if (this.drawMode & DrawMode.SOLID)
		this.model.draw(DrawMode.SOLID);

	// Draw wireframe model?
	if (this.drawMode & DrawMode.WIREFRAME)
		this.model.draw(DrawMode.WIREFRAME);

	// Disable depth test.
	GL.disable(GL.DEPTH_TEST);

	// Draw world axes.
	this.worldAxes.draw(DrawMode.WIREFRAME);

	// Draw model axes.
	this.modelAxes.worldMatrix = this.model.worldMatrix;
	this.modelAxes.draw(DrawMode.WIREFRAME);

	// Draw lights.
	for (var i = 1; i < LIGHT_SOURCE_NUM; i++)
		this.lights[i].draw();
};

// Translate model.
Scene.prototype.translateModel = function(delta) {
	var t = this.camera.volumeToWorld([delta[0], delta[1], delta[2], 0]);
	this.model.translate([t[0], t[1], t[2]]);
};

// Scale model.
Scene.prototype.scaleModel = function(delta) {
	if (delta[1] > 0)
		var s = MODEL_SCALE_UP_FACTOR;
	else
		var s = MODEL_SCALE_DOWN_FACTOR;
	this.model.scale([s, s, s]);
};

// Rotate model.
Scene.prototype.rotateModel = function(delta) {
	// TODO: Why this happenes?
	if (Vec3.norm(delta) == 0)
		return;

	var t = this.camera.volumeToWorld([delta[0], delta[1], delta[2], 0]);

	var n = this.camera.volumeToWorld([0, 0, 1, 0]);
	// TODO: Why this happenes (in perspective mode)?
	if (Vec3.norm(n) == 0)
		n = [0, 0, 1];

	var r = Vec3.cross(Vec3.normalize(t), Vec3.normalize(n));
	var deltaNorm = Vec3.norm(delta);
	this.model.rotate(Vec3.mul(r, [deltaNorm, deltaNorm, deltaNorm]));
};

// Pan camera.
Scene.prototype.panCamera = function(delta) {
	this.camera.pan(delta);
};

// Spin camera.
Scene.prototype.spinCamera = function(delta) {
	this.camera.spin(delta);
};

// Zoom camera.
Scene.prototype.zoomCamera = function(zoomFactor) {
	this.camera.zoom(zoomFactor);
};

/*******************
 *   Class Model   *
 *******************/
 
// Constructor.
Model = function(material) {
	this.visible = true;

	this.objs = new Array();

	this.modelMatrix = Mat4.identity();
	this.worldMatrix = Mat4.identity();
	this.worldTrans = true;
};

// Dump: model

// Initialize model buffers.
Model.prototype.initBuffers = function() {
	for (var i = 0; i < this.objs.length; i++)
		this.objs[i].initBuffers();
};

// Destroy model buffers.
Model.prototype.destroyBuffers = function() {
	for (var i = 0; i < this.objs.length; i++)
		this.objs[i].destroyBuffers();
};

// Draw model.
Model.prototype.draw = function(drawMode) {
	if (!this.visible)
		return;

	// Set matrices.
	this.setMatrices();

	for (var i = 0; i < this.objs.length; i++)
		this.objs[i].draw(drawMode);
};

// Set model matrices.
Model.prototype.setMatrices = function() {
	// Set model matrix.
	var modelMatrix = Mat4.mul(this.modelMatrix, this.worldMatrix);
	GL.uniformMatrix4fv(shaders.loc.uMMat, false, Mat4.toArray(modelMatrix));

	// Set normal matrix.
	var normalMatrix = Mat3.transpose(Mat3.inverse(Mat4.toMat3(modelMatrix)));
	GL.uniformMatrix3fv(shaders.loc.uNMat, false, Mat3.toArray(normalMatrix));
};

// Translate model.
Model.prototype.translate = function(t) {
		if (this.worldTrans)
			this.worldMatrix = Mat4.mul(this.worldMatrix, Mat4.translate(t));
		else
			this.modelMatrix = Mat4.mul(this.modelMatrix, Mat4.translate(t));
};

// Scale model.
Model.prototype.scale = function(s) {
		if (this.worldTrans)
			this.worldMatrix = Mat4.mul(this.worldMatrix, Mat4.scale(s));
		else
			this.modelMatrix = Mat4.mul(this.modelMatrix, Mat4.scale(s));
};

// Rotate model.
Model.prototype.rotate = function(r) {
		if (this.worldTrans)
			this.worldMatrix = Mat4.mul(this.worldMatrix, Mat4.rotateZYX(r));
		else
			this.modelMatrix = Mat4.mul(this.modelMatrix, Mat4.rotateZYX(r));
};

/*************************
 *   Class ModelObject   *
 *************************/

// Constructor.
ModelObject = function(id) {
	this.id = id;

	this.tag = null

	this.RGBA = null;

	this.material = new Object();
	this.material.ambient = MODEL_OBJECT_DEFAULT_MATERIAL_AMBIENT;
	this.material.diffuse = MODEL_OBJECT_DEFAULT_MATERIAL_DIFFUSE;
	this.material.specular = MODEL_OBJECT_DEFAULT_MATERIAL_SPECULAR;
	this.material.emission = MODEL_OBJECT_DEFAULT_MATERIAL_EMISSION;
	this.material.shininess = MODEL_OBJECT_DEFAULT_MATERIAL_SHININESS;

	this.vertices = null;
	this.polys = null;
	this.normals = null;

	this.texture = new Object();
	this.texture.use = false;
	this.texture.index = -1;
	this.texture.UVs = null;
};

// Initialize model object buffers.
ModelObject.prototype.initBuffers = function() {
	// Vertices buffer:
	this.verticesBuffer = GL.createBuffer();
	GL.bindBuffer(GL.ARRAY_BUFFER, this.verticesBuffer);
	GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(this.vertices),
		GL.STATIC_DRAW);

	// Line indices buffer:
	this.lineIndices = [];
	for (var j = 0; j < this.polys.length; j++) {
		var lineIndices = [];
		for (var k = 0; k < this.polys[j].length - 1; k++) {
			lineIndices
				= lineIndices.concat([this.polys[j][k], this.polys[j][k+1]]);
		}
		if (this.tag == ModelObjectTag.POLYGON) {
			lineIndices
				= lineIndices.concat([this.polys[j][k], this.polys[j][0]]);
		}
		this.lineIndices = this.lineIndices.concat(lineIndices);
	}
	this.lineIndicesBuffer = GL.createBuffer();
	GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, this.lineIndicesBuffer);
	GL.bufferData(GL.ELEMENT_ARRAY_BUFFER,
		new Uint16Array(this.lineIndices), GL.STATIC_DRAW);

	// Triangle indices buffer:
	if (this.tag == ModelObjectTag.POLYGON) {
		this.triangleIndices = [];
		for (var j = 0; j < this.polys.length; j++) {
			var triangleIndices = [];
			for (var k = 1; k < this.polys[j].length - 1; k++) {
				triangleIndices = triangleIndices.concat([this.polys[j][0],
					this.polys[j][k], this.polys[j][k+1]]);
			}
			this.triangleIndices
				= this.triangleIndices.concat(triangleIndices);
		}

		this.triangleIndicesBuffer = GL.createBuffer();
		GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, this.triangleIndicesBuffer);
		GL.bufferData(GL.ELEMENT_ARRAY_BUFFER,
			new Uint16Array(this.triangleIndices), GL.STATIC_DRAW);
	}

	// Normals buffer:
	this.normalsBuffer = GL.createBuffer();
	GL.bindBuffer(GL.ARRAY_BUFFER, this.normalsBuffer);
	GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(this.normals),
		GL.STATIC_DRAW);

	// Texture UVs buffer:
	if (this.texture.use) {
		this.texture.UVsBuffer = GL.createBuffer();
		GL.bindBuffer(GL.ARRAY_BUFFER, this.texture.UVsBuffer);
		GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(this.texture.UVs),
			GL.STATIC_DRAW);
	}
};

// Destroy model object buffers.
ModelObject.prototype.destroyBuffers = function() {
	if (this.texture.use)
		destroyBuffer(this.texture.UVsBuffer);
	destroyBuffer(this.normalsBuffer);
	destroyBuffer(this.lineIndicesBuffer);
	if (this.tag == ModelObjectTag.POLYGON)
		destroyBuffer(this.triangleIndicesBuffer);
	destroyBuffer(this.verticesBuffer);
};

// Draw model object.
// TODO: Use display lists when WebGL supports it.
ModelObject.prototype.draw = function(drawMode) {
	// Set model object attributes.
	this.setAttributes(drawMode);

	GL.uniform1i(shaders.loc.uUseLight, 0);
	GL.uniform1i(shaders.loc.uApplyTexture, 0);

	// Draw.
	switch (drawMode) {
	case DrawMode.WIREFRAME:
		GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, this.lineIndicesBuffer);
		GL.drawElements(GL.LINES, this.lineIndices.length,
			GL.UNSIGNED_SHORT, 0);
		break;

	case DrawMode.SOLID:
		if (this.tag == ModelObjectTag.POLYGON) {
			GL.uniform1i(shaders.loc.uUseLight, 1);
			GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, this.triangleIndicesBuffer);
			GL.drawElements(GL.TRIANGLES, this.triangleIndices.length,
				GL.UNSIGNED_SHORT, 0);
		}
		break;

	case DrawMode.TEXTURE:
		if (!this.texture.use)
			return;

		if (this.tag == ModelObjectTag.POLYGON) {
			GL.uniform1i(shaders.loc.uUseLight, 1);

			GL.uniform1i(shaders.loc.uApplyTexture, 1);
			GL.activeTexture(GL.TEXTURE0);
			GL.bindTexture(GL.TEXTURE_2D, textures.data[this.texture.index]);
			GL.uniform1i(shaders.loc.uSampler, 0);

			GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, this.triangleIndicesBuffer);
			GL.drawElements(GL.TRIANGLES, this.triangleIndices.length,
				GL.UNSIGNED_SHORT, 0);
		}
		break;

	case DrawMode.RENDER_TO_TEXTURE:
		GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, this.triangleIndicesBuffer);
		GL.drawElements(GL.TRIANGLES, this.triangleIndices.length,
			GL.UNSIGNED_SHORT, 0);
		break;
	}
};

// Set model object attributes.
ModelObject.prototype.setAttributes = function(drawMode) {
	// Set RGBA.
	var RGBA = this.RGBA;

	if (drawMode == DrawMode.WIREFRAME)
		RGBA = [RGBA[0], RGBA[1], RGBA[2], 1];

	if (drawMode == DrawMode.RENDER_TO_TEXTURE) {
		RGBA = picker.ID2RGBA(this.id);
		for (var i = 0; i < 4; i++)
			RGBA[i] /= 255;
	}

	GL.uniform4fv(shaders.loc.uRGBA, new Float32Array(RGBA));

	var RGB = RGBA.slice(0, 3);

	// Set material.
	var ambient = Vec3.mul(RGB, this.material.ambient);
	GL.uniform3fv(shaders.loc.uMaterial.ambient, new Float32Array(ambient));

	var diffuse = Vec3.mul(RGB, this.material.diffuse);
	GL.uniform3fv(shaders.loc.uMaterial.diffuse, new Float32Array(diffuse));

	var specular = this.material.specular;
	GL.uniform3fv(shaders.loc.uMaterial.specular, new Float32Array(specular));

	var emission = this.material.emission;
	GL.uniform3fv(shaders.loc.uMaterial.emission, new Float32Array(emission));

	var shininess = this.material.shininess;
	GL.uniform1f(shaders.loc.uMaterial.shininess, shininess);

	// Set vertex position.
	GL.bindBuffer(GL.ARRAY_BUFFER, this.verticesBuffer);
	GL.enableVertexAttribArray(shaders.loc.aVertexPosition);
	GL.vertexAttribPointer(shaders.loc.aVertexPosition, 3, GL.FLOAT, false, 0,
		0);

	// Set vertex normal.
	GL.bindBuffer(GL.ARRAY_BUFFER, this.normalsBuffer);
	GL.enableVertexAttribArray(shaders.loc.aVertexNormal);
	GL.vertexAttribPointer(shaders.loc.aVertexNormal, 3, GL.FLOAT, false, 0,
		0);

	// Set vertex texture UV.
	if (this.texture.use) {
		GL.bindBuffer(GL.ARRAY_BUFFER, this.texture.UVsBuffer);
		GL.enableVertexAttribArray(shaders.loc.aVertexTextureUV);
		GL.vertexAttribPointer(shaders.loc.aVertexTextureUV, 2, GL.FLOAT,
			false, 0, 0);
	} else {
		GL.disableVertexAttribArray(shaders.loc.aVertexTextureUV);
	}
};

/********************
 *   Class Camera   *
 ********************/

// Constructor.
Camera = function() {
	// View matrix.
	this.viewAngle = menu.camera.viewAngle;
	this.viewMatrix = Mat4.identity();	// Will be set later ...
	this.viewMatrices = new Array();
	for (var i = 1; i < VIEW_NUM; i++) {
		this.viewMatrices[i] = Mat4.lookat(
			CAMERA_VIEW_EYE[i - 1],
			CAMERA_VIEW_AT[i - 1],
			CAMERA_VIEW_UP[i - 1]
		);
	}

	// Orthographic matrix.
	this.aspectRatio = canvas.width/canvas.height;

	this.left = CAMERA_ORTHO_LEFT * this.aspectRatio;
	this.right = CAMERA_ORTHO_RIGHT * this.aspectRatio;
	this.bottom = CAMERA_ORTHO_BOTTOM;
	this.top = CAMERA_ORTHO_TOP;
	this.near = CAMERA_ORTHO_NEAR;
	this.far = CAMERA_ORTHO_FAR;

	// Orthographic matrix.
	this.orthoMatrix = Mat4.identity();	// Will be set later ...

	// Perspective matrix.
	this.prspMatrix = Mat4.identity();	// Will be set later ...

	// Projection matrix.
	this.projectionMatrix = Mat4.identity();	// Will be set later ...
	
	// Projection mode.
	this.projMode = menu.camera.projMode;
};

// Dump: camera

// Set view.
Camera.prototype.setView = function() {
	this.viewMatrix = this.viewMatrices[this.viewAngle];
};

// Set projection.
Camera.prototype.setProjection = function() {
	this.orthoMatrix = Mat4.ortho(
			this.left, this.right,
			this.bottom, this.top,
			this.near, this.far
		);

	switch (this.projMode) {
	case ProjectionMode.ORTHOGONAL:
		this.projectionMatrix = this.orthoMatrix;
		break;

	case ProjectionMode.PERSPECTIVE:
		this.projectionMatrix = Mat4.mul(this.prspMatrix, this.orthoMatrix);
		break;
	}
};

// Set camera.
Camera.prototype.set = function() {
	this.setMatrices();

	// Set COP (center of projection).
	var COP = Vec4.rMatMul([0, 0, 0, 1], Mat4.inverse(this.viewMatrix)); 
	GL.uniform3fv(shaders.loc.uCOP, new Float32Array(COP.slice(0, 3)));
};

// Set camera matrices.
Camera.prototype.setMatrices = function() {
	// Set view matrix.
	var viewMatrix = this.viewMatrix;
	GL.uniformMatrix4fv(shaders.loc.uVMat, false, Mat4.toArray(viewMatrix));

	// Set projection matrix.
	var projectionMatrix = this.projectionMatrix;
	GL.uniformMatrix4fv(shaders.loc.uPMat, false,
		Mat4.toArray(projectionMatrix));
};

// Volume coordinates to world coordinates.
Camera.prototype.volumeToWorld = function(v) {
	var invViewMat = Mat4.inverse(this.viewMatrix);
	var invProjectionMat = Mat4.inverse(this.projectionMatrix);
	return Vec4.rMatMul(Vec4.rMatMul(v, invProjectionMat), invViewMat);
};

// Pan camera.
Camera.prototype.pan = function(t) {
	this.viewMatrix = Mat4.mul(this.viewMatrix, Mat4.translate(t));
};

// Spin camera.
Camera.prototype.spin = function(r) {
	this.viewMatrix = Mat4.mul(this.viewMatrix, Mat4.rotateZYX(r));
};

// Zoom camera.
Camera.prototype.zoom = function(zoomFactor) {
	this.left *= zoomFactor;
	this.right *= zoomFactor;
	this.bottom *= zoomFactor;
	this.top *= zoomFactor;
	this.setProjection();
};

/*******************
 *   Class Light   *
 *******************/

 // Light type.
 var LightType = {
	AMBIENT: 0,
	POINT: 1,
	DIRECTIONAL: 2,
};

// Constructor.
Light = function(id) {
	this.id = id;

	this.type = null;
	this.TL = null;
	this.pos = null;
	this.dir = null;
	this.attK = null;
	this.ambient = null;
	this.diffuse = null;
	this.specular = null;

	// Lamp:
	this.lamp = new Lamp();

	// Ray:
	this.ray = null;
};

// Set position.
Light.prototype.setPos = function(pos) {
	this.pos = pos;
	this.lamp.worldMatrix
		= Mat4.mul(Mat4.scale([0.1, 0.1, 0.1]), Mat4.translate(this.pos));
};

// Set direction.
Light.prototype.setDir = function(dir) {
	this.dir = dir;

	var nDir = Vec3.normalize(this.dir);
	var nPos = Vec3.cMul(nDir, -1);
	if (this.ray != null)
		this.ray.destroyBuffers();
	this.ray = new Ray(nPos, nDir);
	this.ray.worldMatrix
		= Mat4.mul(Mat4.scale([0.1, 0.1, 0.1]), Mat4.translate(nPos));

	this.ray.visible = this.lamp.visible;
}

// Set light.
Light.prototype.set = function() {
	GL.uniform1i(shaders.loc.uLights[this.id].type, this.type);
	GL.uniform1i(shaders.loc.uLights[this.id].TL, this.TL);

	switch (this.type) {
	case LightType.AMBIENT:
		var ambient = Vec3.cMul(this.ambient, this.ambientI);
		GL.uniform3fv(shaders.loc.uLights[0].ambient,
			new Float32Array(ambient));

		break;

	case LightType.POINT:
		GL.uniform3fv(shaders.loc.uLights[this.id].pos,
			new Float32Array(this.pos));

		GL.uniform3fv(shaders.loc.uLights[this.id].attK,
			new Float32Array(this.attK));

		var diffuse = Vec3.cMul(this.diffuse, this.diffuseI);
		GL.uniform3fv(shaders.loc.uLights[this.id].diffuse,
			new Float32Array(diffuse));

		var specular = Vec3.cMul(this.specular, this.specularI);
		GL.uniform3fv(shaders.loc.uLights[this.id].specular,
			new Float32Array(specular));

		break;

	case LightType.DIRECTIONAL:
		GL.uniform3fv(shaders.loc.uLights[this.id].dir,
			new Float32Array(this.dir));

		var diffuse = Vec3.cMul(this.diffuse, this.diffuseI);
		GL.uniform3fv(shaders.loc.uLights[this.id].diffuse,
			new Float32Array(diffuse));

		var specular = Vec3.cMul(this.specular, this.specularI);
		GL.uniform3fv(shaders.loc.uLights[this.id].specular,
			new Float32Array(specular));

		break;
	}
};

// Draw light.
Light.prototype.draw = function() {
	switch (this.type) {
	case LightType.AMBIENT:
		break;

	case LightType.POINT:
		this.lamp.draw(DrawMode.SOLID);
		break;

	case LightType.DIRECTIONAL:
		this.ray.draw(DrawMode.WIREFRAME);
		break;
	}
};

/******************
 *   Class Axes   *
 ******************/

// Constructor.
Axes = function(RGBA) {
	var axesModel = new Model();

	axesModel.objs[0] = new ModelObject();
	axesModel.objs[0].tag = ModelObjectTag.POLYLINE;
	axesModel.objs[0].RGBA = SCENE_AXES_RGBA;

	axesModel.objs[0].vertices = [
		// x-axis
		0, 0, 0,
		1, 0, 0,
		0.95, 0.01, 0,
		1, 0, 0,
		0.95, -0.01, 0,
		1, 0, 0,
		0.95, -0.03, 0,
		0.97, -0.05, 0, 
		0.95, -0.05, 0,
		0.97, -0.03, 0,

		// y-axis
		0, 0, 0,
		0, 1, 0,
		0.01, 0.95, 0,
		0, 1, 0, 
		-0.01, 0.95, 0,
		0, 1, 0,
		0.03, 0.97, 0,
		0.045, 0.96, 0,
		0.055, 0.97, 0,
		0.035, 0.95, 0,

		// z-axis
		0, 0, 0,
		0, 0, 1,
		0.01, 0, 0.95,
		0, 0, 1,
		-0.01, 0, 0.95,
		0, 0, 1, 
		0.03, 0.01, 1,
		0.05, 0.01, 1,
		0.05, 0.01, 1,
		0.03, -0.01, 1,
		0.03, -0.01, 1,
		0.05, -0.01, 1,
	];

	axesModel.objs[0].normals = [
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
		0, 0, 0,
	];

	axesModel.objs[0].polys = [
		[0, 1],
		[2, 3],
		[4, 5],
		[6, 7],
		[8, 9],

		[10, 11],
		[12, 13],
		[14, 15],
		[16, 17],
		[18, 19],

		[20, 21],
		[22, 23],
		[24, 25],
		[26, 27],
		[28, 29],
		[30, 31],
	];

	axesModel.initBuffers();

	return axesModel;
};

/******************
 *   Class Lamp   *
 ******************/

// Constructor.
Lamp = function() {
	var lampModel = new Model();

	lampModel.visible = false;

	lampModel.objs[0] = new ModelObject();
	lampModel.objs[0].tag = ModelObjectTag.POLYGON;
	lampModel.objs[0].RGBA = LIGHT_LAMP_RGBA;

	lampModel.objs[0].material.ambient = LIGHT_LAMP_MATERIAL_AMBIENT;
	lampModel.objs[0].material.diffuse = LIGHT_LAMP_MATERIAL_DIFFUSE;
	lampModel.objs[0].material.specular = LIGHT_LAMP_MATERIAL_SPECULAR;
	lampModel.objs[0].material.emission = LIGHT_LAMP_MATERIAL_EMISSION;
	lampModel.objs[0].material.shininess = LIGHT_LAMP_MATERIAL_SHININESS;

	lampModel.objs[0].vertices = [
		0.00000, 0.00000, -1.00000,
		-0.38268, 0.00000, -0.92388,
		-0.27060, -0.27060, -0.92388,
		-0.00000, -0.38268, -0.92388,
		0.27060, -0.27060, -0.92388,
		0.38268, 0.00000, -0.92388,
		0.27060, 0.27060, -0.92388,
		-0.00000, 0.38268, -0.92388,
		-0.27060, 0.27060, -0.92388,
		-0.92388, 0.00000, -0.38268,
		-0.65328, -0.65328, -0.38268,
		-0.00000, -0.92388, -0.38268,
		0.65328, -0.65328, -0.38268,
		0.92388, 0.00000, -0.38268,
		0.65328, 0.65328, -0.38268,
		-0.00000, 0.92388, -0.38268,
		-0.65328, 0.65328, -0.38268,
		-0.92388, 0.00000, 0.38268,
		-0.65328, -0.65328, 0.38268,
		-0.00000, -0.92388, 0.38268,
		0.65328, -0.65328, 0.38268,
		0.92388, 0.00000, 0.38268,
		0.65328, 0.65328, 0.38268,
		-0.00000, 0.92388, 0.38268,
		-0.65328, 0.65328, 0.38268,
		-0.38268, 0.00000, 0.92388,
		-0.27060, -0.27060, 0.92388,
		-0.00000, -0.38268, 0.92388,
		0.27060, -0.27060, 0.92388,
		0.38268, 0.00000, 0.92388,
		0.27060, 0.27060, 0.92388,
		-0.00000, 0.38268, 0.92388,
		-0.27060, 0.27060, 0.92388,
		0.38268, -0.00000, 0.92388,
		0.27060, 0.27060, 0.92388,
		0.00000, 0.38268, 0.92388,
		-0.27060, 0.27060, 0.92388,
		-0.38268, 0.00000, 0.92388,
		-0.27060, -0.27060, 0.92388,
		0.00000, -0.38268, 0.92388,
		0.27060, -0.27060, 0.92388,
		0.00000, 0.00000, 1.00000, 
	];

	lampModel.objs[0].normals = [
		0.00000, 0.00000, -1.00000,
		-0.38268, 0.00000, -0.92388,
		-0.27060, -0.27060, -0.92388,
		-0.00000, -0.38268, -0.92388,
		0.27060, -0.27060, -0.92388,
		0.38268, 0.00000, -0.92388,
		0.27060, 0.27060, -0.92388,
		-0.00000, 0.38268, -0.92388,
		-0.27060, 0.27060, -0.92388,
		-0.92388, 0.00000, -0.38268,
		-0.65328, -0.65328, -0.38268,
		-0.00000, -0.92388, -0.38268,
		0.65328, -0.65328, -0.38268,
		0.92388, 0.00000, -0.38268,
		0.65328, 0.65328, -0.38268,
		-0.00000, 0.92388, -0.38268,
		-0.65328, 0.65328, -0.38268,
		-0.92388, 0.00000, 0.38268,
		-0.65328, -0.65328, 0.38268,
		-0.00000, -0.92388, 0.38268,
		0.65328, -0.65328, 0.38268,
		0.92388, 0.00000, 0.38268,
		0.65328, 0.65328, 0.38268,
		-0.00000, 0.92388, 0.38268,
		-0.65328, 0.65328, 0.38268,
		-0.38268, 0.00000, 0.92388,
		-0.27060, -0.27060, 0.92388,
		-0.00000, -0.38268, 0.92388,
		0.27060, -0.27060, 0.92388,
		0.38268, 0.00000, 0.92388,
		0.27060, 0.27060, 0.92388,
		-0.00000, 0.38268, 0.92388,
		-0.27060, 0.27060, 0.92388,
		0.38268, -0.00000, 0.92388,
		0.27060, 0.27060, 0.92388,
		0.00000, 0.38268, 0.92388,
		-0.27060, 0.27060, 0.92388,
		-0.38268, 0.00000, 0.92388,
		-0.27060, -0.27060, 0.92388,
		0.00000, -0.38268, 0.92388,
		0.27060, -0.27060, 0.92388,
		0.00000, 0.00000, 1.00000,
	];

	lampModel.objs[0].polys = [
		[1, 0, 2],
		[2, 0, 3],
		[3, 0, 4],
		[4, 0, 5],
		[5, 0, 6],
		[6, 0, 7],
		[7, 0, 8],
		[8, 0, 1],
		[1, 2, 10, 9],
		[2, 3, 11, 10],
		[3, 4, 12, 11],
		[4, 5, 13, 12],
		[5, 6, 14, 13],
		[6, 7, 15, 14],
		[7, 8, 16, 15],
		[8, 1, 9, 16],
		[9, 10, 18, 17],
		[10, 11, 19, 18],
		[11, 12, 20, 19],
		[12, 13, 21, 20],
		[13, 14, 22, 21],
		[14, 15, 23, 22],
		[15, 16, 24, 23],
		[16, 9, 17, 24],
		[17, 18, 26, 25],
		[18, 19, 27, 26],
		[19, 20, 28, 27],
		[20, 21, 29, 28],
		[21, 22, 30, 29],
		[22, 23, 31, 30],
		[23, 24, 32, 31],
		[24, 17, 25, 32],
		[33, 34, 41],
		[34, 35, 41],
		[35, 36, 41],
		[36, 37, 41],
		[37, 38, 41],
		[38, 39, 41],
		[39, 40, 41],
		[40, 33, 41],
	];

	lampModel.initBuffers();

	return lampModel; 
};

/*****************
 *   Class Ray   *
 *****************/

// Constructor.
Ray = function(from, to) {
	var rayModel = new Model();

	rayModel.visible = false;

	rayModel.objs[0] = new ModelObject();
	rayModel.objs[0].tag = ModelObjectTag.POLYLINE;
	rayModel.objs[0].RGBA = LIGHT_RAY_RGBA;

	rayModel.objs[0].vertices = [
		from[0], from[1], from[2],
		to[0], to[1], to[2],
	];

	rayModel.objs[0].normals = [
		0, 0, 0,
		0, 0, 0,
	];

	rayModel.objs[0].polys = [
		[0, 1],
	];

	rayModel.initBuffers();

	return rayModel; 
};

/******************
 *   Class Vec3   * 
 ******************/

Vec3 = {

	// Clone.
	clone: function(v) {
		return [v[0], v[1], v[2]];
	},

	// Add two 3x1 vectors.
	add: function(u, v) {
		return [u[0] + v[0], u[1] + v[1], u[2] + v[2]];
	},

	// Substract two 3x1 vectors.
	sub: function(u, v) {
		return [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
	},

	// Multiply two 3x1 vectors (element by element).
	mul: function(u, v) {
		return [u[0]*v[0], u[1]*v[1], u[2]*v[2]];
	},

	// Divide two 3x1 vectors (element by element).
	div: function(u, v) {
		return [u[0]/v[0], u[1]/v[1], u[2]/v[2]];
	},

	// Scalar multiplication of two 3x1 vectors.
	dot: function(u, v) {
		return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
	},

	// Cross product of two 3x1 vectors.
	cross: function(u, v) {
		return [
			u[1]*v[2] - u[2]*v[1],
			u[2]*v[0] - u[0]*v[2],
			u[0]*v[1] - u[1]*v[0]
		];
	},

	// Multiply 3x1 vector with a constant.
	cMul: function(v, c) {
		return [c*v[0], c*v[1], c*v[2]];
	},

	// Right multiply 3x1 vector with 3x3 matrix.
	rMatMul: function(v, A) {
		return [
			v[0]*A[0][0] + v[1]*A[1][0] + v[2]*A[2][0],
			v[0]*A[0][1] + v[1]*A[1][1] + v[2]*A[2][1],
			v[0]*A[0][2] + v[1]*A[1][2] + v[2]*A[2][2]
		];
	},

	// Left multiply 3x3 matrix with 3x1 vector.
	lMatMul: function(A, v) {
		return [
			A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2],
			A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2],
			A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2]
		];
	},

	// Norm of a 3x1 vector.
	norm: function(v) {
		return Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	},

	// Normalize the given 3x1 vector.
	normalize: function(v) {
		var norm = this.norm(v);
		return [v[0]/norm, v[1]/norm, v[2]/norm];
	},
};

/******************
 *   Class Vec4   * 
 ******************/

Vec4 = {

	// Clone.
	clone: function(v) {
		return [v[0], v[1], v[2], v[3]];
	},

	// Add two 4x1 vectors.
	add: function(u, v) {
		return [u[0] + v[0], u[1] + v[1], u[2] + v[2], u[3] + v[3]];
	},

	// Substract two 4x1 vectors.
	sub: function(u, v) {
		return [u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3]];
	},

	// Multiply two 4x1 vectors (element by element).
	mul: function(u, v) {
		return [u[0]*v[0], u[1]*v[1], u[2]*v[2], u[3]*v[3]];
	},

	// Divide two 4x1 vectors (element by element).
	div: function(u, v) {
		return [u[0]/v[0], u[1]/v[1], u[2]/v[2], u[3]/v[3]];
	},

	// Multiply 4x1 vector with a constant.
	cMul: function(v, c) {
		return [c*v[0], c*v[1], c*v[2], c*v[3]];
	},

	// Right multiply 4x1 vector with 4x4 matrix.
	rMatMul: function(v, A) {
		return [
			v[0]*A[0][0] + v[1]*A[1][0] + v[2]*A[2][0] + v[3]*A[3][0],
			v[0]*A[0][1] + v[1]*A[1][1] + v[2]*A[2][1] + v[3]*A[3][1],
			v[0]*A[0][2] + v[1]*A[1][2] + v[2]*A[2][2] + v[3]*A[3][2],
			v[0]*A[0][3] + v[1]*A[1][3] + v[2]*A[2][3] + v[3]*A[3][3]
		];
	},

	// Left multiply 4x4 matrix with 4x1 vector.
	lMatMul: function(A, v) {
		return [
			A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2] + A[0][3]*v[3],
			A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2] + A[1][3]*v[3],
			A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2] + A[2][3]*v[3],
			A[3][0]*v[0] + A[3][1]*v[1] + A[3][2]*v[2] + A[3][3]*v[3],
		];
	},

	// Norm of a 4x1 vector.
	norm: function(v) {
		return Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
	},

	// Normalize the given 4x1 vector.
	normalize: function(v) {
		var norm = this.norm(v);
		return [v[0]/norm, v[1]/norm, v[2]/norm, v[3]/norm];
	},
};

/******************
 *   Class Mat3   *
 ******************/

Mat3 = {

	// 3x3 matrix to float 32 array.
	toArray: function(A) {
		return new Float32Array([
			A[0][0], A[0][1], A[0][2],
			A[1][0], A[1][1], A[1][2],
			A[2][0], A[2][1], A[2][2],
		]);
	},

	// 3x3 zero matrix.
	zero: function() {
		return [
			[0, 0, 0],
			[0, 0, 0],
			[0, 0, 0]
		];
	},

	// 3x3 identity matrix.
	identity: function() {
		return [
			[1, 0, 0],
			[0, 1, 0],
			[0, 0, 1]
		];
	},

	// 3x3 random matrix.
	random: function() {
		var A = []
		for (var i = 0; i < 3; i++) {
			A[i] = [];
			for (var j = 0; j < 3; j++)
				A[i][j] = Math.random();
		}
		return A;
	},

	// Add two 3x3 matrices.
	add: function(A, B) {
		var C = [];
		for (var i = 0; i < 3; i++) {
			C[i] = [];
			for (var j = 0; j < 3; j++)
				C[i][j] = A[i][j] + B[i][j];
		}
	},

	// Substract two 3x3 matrices.
	sub: function(A, B) {
		var C = [];
		for (var i = 0; i < 3; i++) {
			C[i] = [];
			for (var j = 0; j < 3; j++)
				C[i][j] = A[i][j] - B[i][j];
		}
	},

	// Multiply two 3x3 matrices.
	mul: function(A, B) {
		var C = this.zero();
		for (var i = 0; i < 3; i++)
			for (var j = 0; j < 3; j++)
				for (var k = 0; k < 3; k++)
					C[i][j] += A[i][k] * B[k][j];
		return C;
	},

	// Multiply 3x3 matrix with a constant.
	cMul: function(A, c) {
		var B = [];
		for (var i = 0; i < 3; i++) {
		B[i] = [];
			for (var j = 0; j < 3; j++)
			B[i][j] = c*A[i][j];
		}
	return B;
	},

	// Transpose 3x3 matrix.
	transpose: function(A) {
		var B = [];
		for (var i = 0; i < 3; i++) {
			B[i] = [];
			for (var j = 0; j < 3; j++)
				B[i][j] = A[j][i];
		}
		return B;
	},

	// Determinant of 3x3 matrix.
	det: function(A) {
		var a11 = A[0][0], a12 = A[0][1], a13 = A[0][2],
			a21 = A[1][0], a22 = A[1][1], a23 = A[1][2],
			a31 = A[2][0], a32 = A[2][1], a33 = A[2][2];

		return a11*a22*a33 - a11*a23*a32
			+ a21*a13*a32 - a21*a12*a33
			+ a31*a12*a23 - a31*a13*a22;
	},

	// Inverse of 3x3 matrix.
	inverse: function(A) {
		var a11 = A[0][0], a12 = A[0][1], a13 = A[0][2],
			a21 = A[1][0], a22 = A[1][1], a23 = A[1][2],
			a31 = A[2][0], a32 = A[2][1], a33 = A[2][2];

		var det = this.det(A);
		var B = new Array();

		B[0] = new Array();
		B[0][0] = (1/det)*(a22*a33 - a23*a32);
		B[0][1] = (1/det)*(a13*a32 - a12*a33);
		B[0][2] = (1/det)*(a12*a23 - a13*a22);

		B[1] = new Array();
		B[1][0] = (1/det)*(a23*a31 - a21*a33);
		B[1][1] = (1/det)*(a11*a33 - a13*a31);
		B[1][2] = (1/det)*(a13*a21 - a11*a23);

		B[2] = new Array();
		B[2][0] = (1/det)*(a21*a32 - a22*a31);
		B[2][1] = (1/det)*(a12*a31 - a11*a32);
		B[2][2] = (1/det)*(a11*a22 - a12*a21);

		return B;
	}
};

/******************
 *   Class Mat4   *
 ******************/

Mat4 = {

	// 4x4 matrix to float 32 array.
	toArray: function(A) {
		return new Float32Array([
			A[0][0], A[0][1], A[0][2], A[0][3],
			A[1][0], A[1][1], A[1][2], A[1][3],
			A[2][0], A[2][1], A[2][2], A[2][3],
			A[3][0], A[3][1], A[3][2], A[3][3],
		]);
	},

	// 3x3 zero matrix.
	zero: function() {
		return [
			[0, 0, 0, 0],
			[0, 0, 0, 0],
			[0, 0, 0, 0],
			[0, 0, 0, 0]
		];
	},

	// 3x3 identity matrix.
	identity: function() {
		return [
			[1, 0, 0, 0],
			[0, 1, 0, 0],
			[0, 0, 1, 0],
			[0, 0, 0, 1]
		];
	},

	// 4x4 random matrix.
	random: function() {
		var A = []
		for (var i = 0; i < 4; i++) {
			A[i] = [];
			for (var j = 0; j < 4; j++)
				A[i][j] = Math.random();
		}
		return A;
	},

	// Add two 4x4 matrices.
	add: function(A, B) {
		var C = [];
		for (var i = 0; i < 4; i++) {
			C[i] = [];
			for (var j = 0; j < 4; j++)
				C[i][j] = A[i][j] + B[i][j];
		}
	},

	// Substract two 4x4 matrices.
	sub: function(A, B) {
		var C = [];
		for (var i = 0; i < 4; i++) {
			C[i] = [];
			for (var j = 0; j < 4; j++)
				C[i][j] = A[i][j] - B[i][j];
		}
	},

	// Multiply two 4x4 matrices.
	mul: function(A, B) {
		var C = this.zero();
		for (var i = 0; i < 4; i++)
			for (var j = 0; j < 4; j++)
				for (var k = 0; k < 4; k++)
					C[i][j] += A[i][k] * B[k][j];
		return C;
	},

	// Multiply 4x4 matrix with a constant.
	cMul: function(A, c) {
		var B = [];
		for (var i = 0; i < 4; i++) {
			B[i] = [];
			for (var j = 0; j < 4; j++)
				B[i][j] = c*A[i][j];
		}
		return B;
	},

	// Transpose 4x4 matrix.
	transpose: function(A) {
		var B = [];
		for (var i = 0; i < 4; i++) {
			B[i] = [];
			for (var j = 0; j < 4; j++)
				B[i][j] = A[j][i];
		}
		return B;
	},

	// Determinant of 4x4 matrix.
	det: function(A) {
		var a11 = A[0][0], a12 = A[0][1], a13 = A[0][2], a14 = A[0][3],
			a21 = A[1][0], a22 = A[1][1], a23 = A[1][2], a24 = A[1][3],
			a31 = A[2][0], a32 = A[2][1], a33 = A[2][2], a34 = A[2][3],
			a41 = A[3][0], a42 = A[3][1], a43 = A[3][2], a44 = A[3][3];
	
		return a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43
			+ a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41
			+ a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42
			+ a14*a21*a33*a42 + a14*a22*a31*a43 + a14*a23*a32*a41
			- a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42
			- a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43
			- a13*a21*a34*a42 - a13*a22*a31*a44 - a13*a24*a32*a41
			- a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42;
	},

	// Inverse of 4x4 matrix.
	inverse: function(A) {
		var a11 = A[0][0], a12 = A[0][1], a13 = A[0][2], a14 = A[0][3],
			a21 = A[1][0], a22 = A[1][1], a23 = A[1][2], a24 = A[1][3],
			a31 = A[2][0], a32 = A[2][1], a33 = A[2][2], a34 = A[2][3],
			a41 = A[3][0], a42 = A[3][1], a43 = A[3][2], a44 = A[3][3];
	
		var det = this.det(A);
		var B = new Array();
	
		B[0] = new Array();
		B[0][0] = (1/det)*(a22*a33*a44 + a23*a34*a42 + a24*a32*a43
			- a22*a34*a43 - a23*a32*a44 - a24*a33*a42);
		B[0][1] = (1/det)*(a12*a34*a43 + a13*a32*a44 + a14*a33*a42
			- a12*a33*a44 - a13*a34*a42 - a14*a32*a43);
		B[0][2] = (1/det)*(a12*a23*a44 + a13*a24*a42 + a14*a22*a43
			- a12*a24*a43 - a13*a22*a44 - a14*a23*a42);
		B[0][3] = (1/det)*(a12*a24*a33 + a13*a22*a34 + a14*a23*a32
			- a12*a23*a34 - a13*a24*a32 - a14*a22*a33);
		
		B[1] = new Array();
		B[1][0] = (1/det)*(a21*a34*a43 + a23*a31*a44 + a24*a33*a41
			- a21*a33*a44 - a23*a34*a41 - a24*a31*a43);
		B[1][1] = (1/det)*(a11*a33*a44 + a13*a34*a41 + a14*a31*a43
			- a11*a34*a43 - a13*a31*a44 - a14*a33*a41);
		B[1][2] = (1/det)*(a11*a24*a43 + a13*a21*a44 + a14*a23*a41
			- a11*a23*a44 - a13*a24*a41 - a14*a21*a43);
		B[1][3] = (1/det)*(a11*a23*a34 + a13*a24*a31 + a14*a21*a33
			- a11*a24*a33 - a13*a21*a34 - a14*a23*a31);
		
		B[2] = new Array();
		B[2][0] = (1/det)*(a21*a32*a44 + a22*a34*a41 + a24*a31*a42
			- a21*a34*a42 - a22*a31*a44 - a24*a32*a41);
		B[2][1] = (1/det)*(a11*a34*a42 + a12*a31*a44 + a14*a32*a41
			- a11*a32*a44 - a12*a34*a41 - a14*a31*a42);
		B[2][2] = (1/det)*(a11*a22*a44 + a12*a24*a41 + a14*a21*a42
			- a11*a24*a42 - a12*a21*a44 - a14*a22*a41);
		B[2][3] = (1/det)*(a11*a24*a32 + a12*a21*a34 + a14*a22*a31
			- a11*a22*a34 - a12*a24*a31 - a14*a21*a32);
		
		B[3] = new Array();
		B[3][0] = (1/det)*(a21*a33*a42 + a22*a31*a43 + a23*a32*a41
			- a21*a32*a43 - a22*a33*a41 - a23*a31*a42);
		B[3][1] = (1/det)*(a11*a32*a43 + a12*a33*a41 + a13*a31*a42
			-a11*a33*a42 - a12*a31*a43 - a13*a32*a41);
		B[3][2] = (1/det)*(a11*a23*a42 + a12*a21*a43 + a13*a22*a41
			- a11*a22*a43 - a12*a23*a41 - a13*a21*a42);
		B[3][3] = (1/det)*(a11*a22*a33 + a12*a23*a31 + a13*a21*a32
			- a11*a23*a32 - a12*a21*a33 - a13*a22*a31);
	
		return B;
	},

	// Translation matrix.
	translate: function(t) {
		return [
			[   1,    0,    0, 0],
			[   0,    1,    0, 0],
			[   0,    0,    1, 0],
			[t[0], t[1], t[2], 1]
		];
	},

	// Scale matrix.
	scale: function(s) {
		return [
			[s[0],    0,    0, 0],
			[   0, s[1],    0, 0],
			[   0,    0, s[2], 0],
			[   0,    0,    0, 1]
		];
	},

	// Rotation around x-axis matrix.
	rotateX: function(theta) {
		var ct = Math.cos(theta);
		var st = Math.sin(theta);
		return [
			[1,   0,  0, 0],
			[0,  ct, st, 0],
			[0, -st, ct, 0],
			[0,   0,  0, 1]
		];
	},

	// Rotation around y-axis matrix.
	rotateY: function(theta) {
		var ct = Math.cos(theta);
		var st = Math.sin(theta);
		return [
			[ct, 0, -st, 0],
			[ 0, 1,   0, 0],
			[st, 0,  ct, 0],
			[ 0, 0,   0, 1]
		];
	},

	// Rotation around z-axis matrix.
	rotateZ: function(theta) {
		var ct = Math.cos(theta);
		var st = Math.sin(theta);
		return [
			[ ct, st, 0, 0],
			[-st, ct, 0, 0],
			[  0,  0, 1, 0],
			[  0,  0, 0, 1]
		];
	},

	// Rotation around x,y,z axes matrix.
	rotateZYX: function(r) {
		return this.mul(
			this.mul(this.rotateZ(r[2]), this.rotateY(r[1])),
			this.rotateX(r[0])
		);
	},

	// Look at matrix.
	lookat: function(eye, at, up) {
		var n = Vec3.normalize(Vec3.sub(eye, at));
		var u = Vec3.normalize(Vec3.cross(up, n));
		var v = Vec3.cross(n, u);
	
		var C = [
			[u[0], u[1], u[2], 0],
			[v[0], v[1], v[2], 0],
			[n[0], n[1], n[2], 0],
			[0, 0, 0, 1]
		];

		var T = this.translate(Vec3.cMul(eye, -1));
	
		return this.mul(T, this.transpose(C));
	},

	// Orthographic matrix.
	ortho: function(left, right, bottom, top, near, far) {
		O = this.zero();	
		O[0][0] = 2 / (right - left);
		O[3][0] = -(right + left) / (right - left);
		O[1][1] = 2 / (top - bottom);
		O[3][1] = -(top + bottom) / (top - bottom);
		O[2][2] = -2 / (far - near);
		O[3][2] = -(far + near) / (far - near);
		O[3][3] = 1;

		return O;
	},

	// Frustum matrix.
	frustum: function(left, right, bottom, top, near, far) {
		F = this.zero();
		F[0][0] = 2 * near / (right - left);
		F[2][0] = (right + left) / (right - left);
		F[1][1] = 2 * near / (top - bottom);
		F[2][1] = (top + bottom) / (top - bottom);
		F[2][2] = -(far + near) / (far - near);
		F[3][2] = -2 * far * near / (far - near);
		F[2][3] = -1;

		return F;
	},

	// Perspective matrix.
	perspective: function(fovy, aspect, near, far) { 
		var top = near * Math.tan(fovy * Math.PI / 180);
		var right = top * aspect;
		var left = -right;
		var bottom = -top; 
		return this.frustum(left, right, bottom, top, near, far);
	},

	// 4x4 matrix to 3x3 matrix.
	toMat3: function(A) {
		return [
		[A[0][0], A[0][1], A[0][2]],
		[A[1][0], A[1][1], A[1][2]],
		[A[2][0], A[2][1], A[2][2]]
		];
	}
};
