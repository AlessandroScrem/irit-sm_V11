/*****************************************************************************
* A filter to convert Viewer JS encoded string to a JS file.		     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber/Jesse Oberstein		Ver 2.0, July 2014   *
*****************************************************************************/

#include "irit23js.h"

IRIT_GLOBAL_DATA char
    *Irit23jsGlblViewerString = "\
    <!DOCTYPE html>\n\
    <html lang=\"en\">\n\
    <head>\n\
    <meta charset=\"utf-8\">\n\
	<title>Irit to Three.js</title>\n\
    </head>\n\
    <body style='font-family: Verdana'>\n\
	<div class='reset' style='float: right; border: 2px solid #000; padding: 5px 10px 5px;'>\n\
	    <p>Reset Options:</p>\n\
	    <button type=\"button\" onclick=\"resetCamera()\">\n\
		Reset Camera</button>\n\
	    <button type=\"button\" onclick=\"culling = 1, wires = false, determineRenderSettings()\">\n\
		Default Settings</button>\n\
	</div>\n\
	\n\
	<div class='culling' style='float: right; border: 2px solid #000; margin-right: -2px; padding: 5px 10px 5px;'>\n\
	    <p>Culling Options:</p>\n\
	    <button type=\"button\" onclick=\"culling = 0, determineRenderSettings()\">\n\
		    Backface Culling</button>\n\
	    <button type=\"button\" onclick=\"culling = 1, determineRenderSettings()\">\n\
		    No Culling</button>\n\
	    <button type=\"button\" onclick=\"culling = 2, determineRenderSettings()\">\n\
		    Frontface Culling</button>\n\
	</div>\n\
	\n\
	<div class='wireframe' style='float: right; border: 2px solid #000; margin-right: -2px; padding: 5px 10px 5px;'>\n\
	    <p>Wireframe Options:</p>\n\
	    <button type=\"button\" onclick=\"wires = true, determineRenderSettings()\">\n\
		    Turn On Wireframe</button>\n\
	    <button type=\"button\" onclick=\"wires = false, determineRenderSettings()\">\n\
		    Turn Off Wireframe</button>\n\
	</div>\n\
	<h1 style='margin-left: 20px'> Irit to THREE.js Filter</h1>\n\
	<h4 style='margin-left: 60px; font-size: 15px;'> Camera: %s || File converted: %s</h4>\n\
	<div id='viewer-back' style='margin: -10px auto 0; background-color: #B0E0E6; border: 1px solid #0000FF;'>\n\
	    <div id='viewer' style='margin: 0 auto 0;'></div>\n\
	</div>\n\
	<script type='text/javascript' src=\"irit3js.js\"></script>\n\
	<script type='text/javascript' src=\"iritOC.js\"></script>\n\
	<script type='text/javascript'>\n\
		// Variables used to set up the scene.\n\
		var camera, scene, ambientLight, directionalLight, renderer, controls;\n\
		// The default culling option (culling is disabled).\n\
		var culling = 1;\n\
		// The material to be placed on the loaded object. \n\
		var material; \n\
		// The default wireframe option (wireframe is disabled). \n\
		var wires = false; \n\
		// The rate to zoom in or out at.  Used only with the Orthographic camera.\n\
		var zoom = 1;\n\
		\n\
		// Sets a new scene for the viewer.\n\
		function displayScene() {\n\
		    if (scene == undefined) {\n\
			    scene = new THREE.Scene();\n\
			    %s\n\
			    camera.position.x = 0;\n\
			    camera.position.y = 0;\n\
			    camera.position.z = 7;\n\
			    \n\
			    ambientLight = new THREE.AmbientLight(0x555555);\n\
			    scene.add(ambientLight);\n\
			    \n\
			    directionalLight = new THREE.DirectionalLight(0xffffff);\n\
			    directionalLight.position.set(1, 1, 1).normalize();\n\
			    scene.add(directionalLight);\n\
			    \n\
			    renderer = new THREE.WebGLRenderer();\n\
			    var w = window.innerWidth\n\
				|| document.documentElement.clientWidth\n\
				|| document.body.clientWidth;\n\
			    \n\
			    var h = window.innerHeight\n\
				|| document.documentElement.clientHeight\n\
				|| document.body.clientHeight;\n\
				w *= .86;\n\
				h *= .86;\n\
			    renderer.setSize(w, h);\n\
			    renderer.setClearColor(0xB0E0E6);\n\
			    %s\n\
			    \n\
			    document.getElementById('viewer').style.width = w.toString().concat('px');\n\
			    document.getElementById('viewer').style.height = h.toString().concat('px');\n\
			    document.getElementById('viewer').appendChild(renderer.domElement);\n\
			    controls = new THREE.OrbitControls(camera, renderer.domElement);\n\
		    }\n\
		    else { \n\
			    scene = new THREE.Scene(); \n\
			    scene.add(ambientLight); \n\
			    scene.add(directionalLight); \n\
		    }\n\
		}\n\
		\n\
		// Allows a user to manipulate the view.\n\
		function animate() {\n\
		    requestAnimationFrame(animate);\n\
		    renderer.render(scene, camera);\n\
		    controls.update();\n\
		}\n\
		\n\
		%s\n\
		\n\
		// Resets the camera and model to their original positions.\n\
		function resetCamera() {\n\
			%s\n\
		}\n\
		\n\
		// Determines the culling, sets the wrapping if a texture exists, and sets the shading.\n\
		function determineRenderSettings() {\n\
		    for (var i = 0; i < material.materials.length; i++) {\n\
			    // Enables backface culling. \n\
			    if (culling == 0) {\n\
				material.materials[i].side = THREE.FrontSide;\n\
			    }\n\
			    // Disables culling. \n\
			    else if (culling == 1) {\n\
				material.materials[i].side = THREE.DoubleSide;\n\
			    }\n\
			    // Enables frontface culling. \n\
			    else if (culling == 2) {\n\
				material.materials[i].side = THREE.BackSide;\n\
			    }\n\
			    // Sets a texture map is a texture exists for this model. \n\
			    if (material.materials[i].map != undefined) { \n\
				%s \n\
				material.materials[i].map.wrapS = THREE.RepeatWrapping; \n\
				material.materials[i].map.wrapT = THREE.RepeatWrapping; \n\
			    } \n\
			    // Enables wireframe mode. \n\
			    if (wires) { \n\
				material.materials[i].wireframe = true;\n\
			    }\n\
			    // Disables frontface culling. \n\
			    if (!wires) { \n\
				material.materials[i].wireframe = false;\n\
			    }\n\
			    material.materials[i].shading = THREE.SmoothShading;\n\
		    }\n\
		}\n\
		\n\
		// Loads the provided JSON file.\n\
		function displayJSONFile() {\n\
		    displayScene();\n\
		    var loader = new THREE.JSONLoader();\n\
		    loader.load(\"%s\", function (geometry, mtrl) {\n\
			    material = new THREE.MeshFaceMaterial( mtrl );\n\
			    var mesh = new THREE.Mesh(geometry, material);\n\
			    determineRenderSettings(); \n\
			    scene.add(mesh);\n\
		    });\n\
		    animate();\n\
		}\n\
		\n\
		// Does the user's browser support WebGL?\n\
		function supportsWebGL() {\n\
			try { \n\
				var canvas = document.createElement( 'canvas' ); \n\
				return !! window.WebGLRenderingContext && ( canvas.getContext( 'webgl' ) \n\
					|| canvas.getContext( 'experimental-webgl' ) ); \n\
			} \n\
			catch( e ) { \n\
				return false; \n\
			} \n\
		}\n\
		\n\
		// Calls the window to display the JSON file if WebGL is supported.\n\
		if ( ! supportsWebGL() ) {\n\
		    alert(\"Sorry, WebGL is not supported by this browser. Please upgrade to the newest version of your current browser, or use Google Chrome or Firefox.\");\n\
		}\n\
		else {\n\
		    displayJSONFile();\n\
		}\n\
		\n\
	</script>\n\
    </body>\n\
    </html>";
