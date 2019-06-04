/*****************************************************************************
* A filter to convert Orbit control JS encoded string to a JS file.	     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber/Jesse Oberstein		Ver 2.0, July 2014   *
*****************************************************************************/

#include "irit23js.h"

/* The following JS code is from: http://threejs.org/ (downlaod) */

IRIT_GLOBAL_DATA char
    *Irit23jsGlblOCString[] = { "\
	THREE.OrbitControls = function ( object, domElement ) {\n\
	this.object = object;\n\
	this.domElement = ( domElement !== undefined ) ? domElement : document;\n\
	\n\
	// API\n\
	\n\
	// Set to false to disable this control\n\
	this.enabled = true;\n\
	\n\
	// \"target\" sets the location of focus, where the control orbits around\n\
	// and where it pans with respect to.\n\
	this.target = new THREE.Vector3();\n\
	\n\
	// center is old, deprecated; use \"target\" instead\n\
	this.center = this.target;\n\
	\n\
	// This option actually enables dollying in and out; left as \"zoom\" for\n\
	// backwards compatibility\n\
	this.noZoom = false;\n\
	this.zoomSpeed = 1.0;\n\
	\n\
	// Limits to how far you can dolly in and out\n\
	this.minDistance = 0;\n\
	this.maxDistance = Infinity;\n\
	\n\
	// Set to true to disable this control\n\
	this.noRotate = false;\n\
	this.rotateSpeed = 1.0;\n\
	\n\
	// Set to true to disable this control\n\
	this.noPan = false;\n\
	this.keyPanSpeed = 7.0;	// pixels moved per arrow key push\n\
	\n\
	// Set to true to automatically rotate around the target\n\
	this.autoRotate = false;\n\
	this.autoRotateSpeed = 2.0; // 30 seconds per round when fps is 60\n\
	\n\
	// How far you can orbit vertically, upper and lower limits.\n\
	// Range is 0 to Math.PI radians.\n\
	this.minPolarAngle = 0; // radians\n\
	this.maxPolarAngle = Math.PI; // radians\n\
	\n\
	// Set to true to disable use of the keys\n\
	this.noKeys = false;\n\
	\n\
	// The four arrow keys\n\
	this.keys = { LEFT: 37, UP: 38, RIGHT: 39, BOTTOM: 40 };\n\
	\n\
	////////////\n\
	// internals\n\
	\n\
	var scope = this;\n\
	\n\
	var EPS = 0.000001;\n\
	\n\
	var rotateStart = new THREE.Vector2();\n\
	var rotateEnd = new THREE.Vector2();\n\
	var rotateDelta = new THREE.Vector2();\n\
	\n\
	var panStart = new THREE.Vector2();\n\
	var panEnd = new THREE.Vector2();\n\
	var panDelta = new THREE.Vector2();\n\
	var panOffset = new THREE.Vector3();\n\
	\n\
	var offset = new THREE.Vector3();\n\
	\n\
	var dollyStart = new THREE.Vector2();\n\
	var dollyEnd = new THREE.Vector2();\n\
	var dollyDelta = new THREE.Vector2();\n\
	\n\
	var phiDelta = 0;\n\
	var thetaDelta = 0;\n\
	var scale = 1;\n\
	var pan = new THREE.Vector3();\n\
	\n\
	var lastPosition = new THREE.Vector3();\n\
	\n\
	var STATE = { NONE : -1, ROTATE : 0, DOLLY : 1, PAN : 2, TOUCH_ROTATE : 3, TOUCH_DOLLY : 4, TOUCH_PAN : 5 };\n\
	\n\
	var state = STATE.NONE;\n\
	\n\
	// for reset\n\
	\n\
	this.target0 = this.target.clone();\n\
	this.position0 = this.object.position.clone();\n\
	\n\
	// so camera.up is the orbit axis\n\
	\n\
	var quat = new THREE.Quaternion().setFromUnitVectors( object.up, new THREE.Vector3( 0, 1, 0 ) );\n\
	var quatInverse = quat.clone().inverse();\n\
	\n\
	// events\n\
	\n\
	var changeEvent = { type: 'change' };\n\
	var startEvent = { type: 'start'};\n\
	var endEvent = { type: 'end'};\n\
	\n\
	this.rotateLeft = function ( angle ) {\n\
	\n\
		if ( angle === undefined ) {\n\
		\n\
			angle = getAutoRotationAngle();\n\
			\n\
		}\n\
		\n\
		thetaDelta -= angle;\n\
		\n\
	};\n\
	\n\
	this.rotateUp = function ( angle ) {\n\
	\n\
		if ( angle === undefined ) {\n\
		\n\
			angle = getAutoRotationAngle();\n\
			\n\
		}\n\
		\n\
		phiDelta -= angle;\n\
		\n\
	};\n\
	\n\
	// pass in distance in world space to move left\n\
	this.panLeft = function ( distance ) {\n\
	\n\
		var te = this.object.matrix.elements;\n\
		\n\
		// get X column of matrix\n\
		panOffset.set( te[ 0 ], te[ 1 ], te[ 2 ] );\n\
		panOffset.multiplyScalar( - distance );\n\
		\n\
		pan.add( panOffset );\n\
		\n\
	};\n\
	\n\
	// pass in distance in world space to move up\n\
	this.panUp = function ( distance ) {\n\
	\n\
		var te = this.object.matrix.elements;\n\
		\n\
		// get Y column of matrix\n\
		panOffset.set( te[ 4 ], te[ 5 ], te[ 6 ] );\n\
		panOffset.multiplyScalar( distance );\n\
		\n\
		pan.add( panOffset );\n\
		\n\
	};\n\
	\n\
	// pass in x,y of change desired in pixel space,\n\
	// right and down are positive\n\
	this.pan = function ( deltaX, deltaY ) {\n\
	\n\
		var element = scope.domElement === document ? scope.domElement.body : scope.domElement;\n\
		\n\
		if ( scope.object.fov !== undefined ) {\n\
		\n\
			// perspective\n\
			var position = scope.object.position;\n\
			var offset = position.clone().sub( scope.target );\n\
			var targetDistance = offset.length();\n\
			\n\
			// half of the fov is center to top of screen\n\
			targetDistance *= Math.tan( ( scope.object.fov / 2 ) * Math.PI / 180.0 );\n\
			\n\
			// we actually don't use screenWidth, since perspective camera is fixed to screen height\n\
			scope.panLeft( 2 * deltaX * targetDistance / element.clientHeight );\n\
			scope.panUp( 2 * deltaY * targetDistance / element.clientHeight );\n\
			\n\
		} else if ( scope.object.top !== undefined ) {\n\
		\n\
			// orthographic\n\
			scope.panLeft( deltaX * (scope.object.right - scope.object.left) / element.clientWidth );\n\
			scope.panUp( deltaY * (scope.object.top - scope.object.bottom) / element.clientHeight );\n\
			\n\
		} else {\n\
		\n\
			// camera neither orthographic or perspective\n\
			console.warn( 'WARNING: OrbitControls.js encountered an unknown camera type - pan disabled.' );\n\
			\n\
		}\n\
		\n\
	};\n\
	\n\
	this.dollyIn = function ( dollyScale ) {\n\
	\n\
		if ( dollyScale === undefined ) {\n\
		\n\
			dollyScale = getZoomScale();\n\
			\n\
		}\n\
		\n\
		scale /= dollyScale;\n\
		\n\
	};\n\
	\n\
	this.dollyOut = function ( dollyScale ) {\n\
	\n\
		if ( dollyScale === undefined ) {\n\
		\n\
			dollyScale = getZoomScale();\n\
			\n\
		}\n\
		\n\
		scale *= dollyScale;\n\
		\n\
	};\n\
	\n\
	this.update = function () {\n\
	\n\
		var position = this.object.position;\n\
		\n\
		offset.copy( position ).sub( this.target );\n\
		\n\
		// rotate offset to \"y-axis-is-up\" space\n\
		offset.applyQuaternion( quat );\n\
		\n\
		// angle from z-axis around y-axis\n\
		\n\
		var theta = Math.atan2( offset.x, offset.z );\n\
		\n\
		// angle from y-axis\n\
		\n\
		var phi = Math.atan2( Math.sqrt( offset.x * offset.x + offset.z * offset.z ), offset.y );\n\
		\n\
		if ( this.autoRotate ) {\n\
		\n\
			this.rotateLeft( getAutoRotationAngle() );\n\
			\n\
		}\n\
		\n\
		theta += thetaDelta;\n\
		phi += phiDelta;\n\
		\n\
		// restrict phi to be between desired limits\n\
		phi = Math.max( this.minPolarAngle, Math.min( this.maxPolarAngle, phi ) );\n\
		\n\
		// restrict phi to be betwee EPS and PI-EPS\n\
		phi = Math.max( EPS, Math.min( Math.PI - EPS, phi ) );\n\
		\n\
		var radius = offset.length() * scale;\n\
		\n\
		// restrict radius to be between desired limits\n\
		radius = Math.max( this.minDistance, Math.min( this.maxDistance, radius ) );\n\
		\n\
		// move target to panned location\n\
		this.target.add( pan );\n\
		\n\
		offset.x = radius * Math.sin( phi ) * Math.sin( theta );\n\
		offset.y = radius * Math.cos( phi );\n\
		offset.z = radius * Math.sin( phi ) * Math.cos( theta );\n\
		\n\
		// rotate offset back to \"camera-up-vector-is-up\" space\n\
		offset.applyQuaternion( quatInverse );\n\
		\n\
		position.copy( this.target ).add( offset );\n\
		\n\
		this.object.lookAt( this.target );\n\
		\n\
		thetaDelta = 0;\n\
		phiDelta = 0;\n\
		scale = 1;\n\
		pan.set( 0, 0, 0 );\n\
		\n\
		if ( lastPosition.distanceToSquared( this.object.position ) > EPS ) {\n\
		\n\
			this.dispatchEvent( changeEvent );\n\
			\n\
			lastPosition.copy( this.object.position );\n\
			\n\
		}\n\
		\n\
	};\n\
	\n\
	this.reset = function () {\n\
	\n\
		state = STATE.NONE;\n\
		\n\
		this.target.copy( this.target0 );\n\
		this.object.position.copy( this.position0 );\n\
		\n\
		this.update();\n\
		\n\
	};\n\
	\n\
	function getAutoRotationAngle() {\n\
	\n\
		return 2 * Math.PI / 60 / 60 * scope.autoRotateSpeed;\n\
		\n\
	}\n\
	\n\
	function getZoomScale() {\n\
	\n\
		return Math.pow( 0.95, scope.zoomSpeed );\n\
		\n\
	}\n\
	\n\
	function onMouseDown( event ) {\n\
	\n\
		if ( scope.enabled === false ) return;\n\
		event.preventDefault();\n\
		\n\
		if ( event.button === 0 ) {\n\
			if ( scope.noRotate === true ) return;\n\
			\n\
			state = STATE.ROTATE;\n\
			\n\
			rotateStart.set( event.clientX, event.clientY );\n\
			\n\
		} else if ( event.button === 1 ) {\n\
			if ( scope.noZoom === true ) return;\n\
			\n\
			state = STATE.DOLLY;\n\
			\n\
			dollyStart.set( event.clientX, event.clientY );\n\
			\n\
		} else if ( event.button === 2 ) {\n\
			if ( scope.noPan === true ) return;\n\
			\n\
			state = STATE.PAN;\n\
			\n\
			panStart.set( event.clientX, event.clientY );\n\
			\n\
		}\n\
		\n\
		scope.domElement.addEventListener( 'mousemove', onMouseMove, false );\n\
		scope.domElement.addEventListener( 'mouseup', onMouseUp, false );\n\
		scope.dispatchEvent( startEvent );\n\
		\n\
	}\n\
	\n\
	function onMouseMove( event ) {\n\
	\n\
		if ( scope.enabled === false ) return;\n\
		\n\
		event.preventDefault();\n\
		\n\
		var element = scope.domElement === document ? scope.domElement.body : scope.domElement;\n\
		\n\
		if ( state === STATE.ROTATE ) {\n\
		\n\
			if ( scope.noRotate === true ) return;\n\
			\n\
			rotateEnd.set( event.clientX, event.clientY );\n\
			rotateDelta.subVectors( rotateEnd, rotateStart );\n\
			\n\
			// rotating across whole screen goes 360 degrees around\n\
			scope.rotateLeft( 2 * Math.PI * rotateDelta.x / element.clientWidth * scope.rotateSpeed );\n\
			\n\
			// rotating up and down along whole screen attempts to go 360, but limited to 180\n\
			scope.rotateUp( 2 * Math.PI * rotateDelta.y / element.clientHeight * scope.rotateSpeed );\n\
			\n\
			rotateStart.copy( rotateEnd );\n\
			\n\
		} else if ( state === STATE.DOLLY ) {\n\
		\n\
			if ( scope.noZoom === true ) return;\n\
			\n\
			dollyEnd.set( event.clientX, event.clientY );\n\
			dollyDelta.subVectors( dollyEnd, dollyStart );\n\
			\n\
			if ( dollyDelta.y > 0 ) {\n\
			\n\
				scope.dollyIn();\n\
				\n\
			} else {\n\
			\n\
				scope.dollyOut();\n\
				\n\
			}\n\
			dollyStart.copy( dollyEnd );\n\
			\n\
		} else if ( state === STATE.PAN ) {\n\
		\n\
			if ( scope.noPan === true ) return;\n\
			\n\
			panEnd.set( event.clientX, event.clientY );\n\
			panDelta.subVectors( panEnd, panStart );\n\
			\n\
			scope.pan( panDelta.x, panDelta.y );\n\
			\n\
			panStart.copy( panEnd );\n\
			\n\
		}\n\
		\n\
		scope.update();\n\
		\n\
	}\n\
	\n\
	function onMouseUp( ) { \n\
	\n\
		if ( scope.enabled === false ) return;\n\
		\n\
		scope.domElement.removeEventListener( 'mousemove', onMouseMove, false );\n\
		scope.domElement.removeEventListener( 'mouseup', onMouseUp, false );\n\
		scope.dispatchEvent( endEvent );\n\
		state = STATE.NONE;\n\
		\n\
	}\n\
	\n\
	function onMouseWheel( event ) {\n\
	\n\
		if ( scope.enabled === false || scope.noZoom === true ) return;\n\
		\n\
		event.preventDefault();\n\
		event.stopPropagation();\n\
		\n\
		var delta = 0;\n\
		\n\
		if ( event.wheelDelta !== undefined ) { // WebKit / Opera / Explorer 9\n\
		\n\
			delta = event.wheelDelta;\n\
			\n\
		} else if ( event.detail !== undefined ) { // Firefox\n\
		\n\
			delta = - event.detail;\n\
			\n\
		}\n\
		\n\
		if ( delta > 0 ) {\n\
		\n\
			scope.dollyOut();\n\
			\n\
		} else {\n\
		\n\
			scope.dollyIn();\n\
			\n\
		}\n\
		\n\
		scope.update();\n\
		scope.dispatchEvent( startEvent );\n\
		scope.dispatchEvent( endEvent );\n\
		\n\
	}\n\
	\n\
	function onKeyDown( event ) {\n\
	\n\
		if ( scope.enabled === false || scope.noKeys === true || scope.noPan === true ) return;\n\
		\n\
		switch ( event.keyCode ) {\n\
		\n\
			case scope.keys.UP:\n\
				scope.pan( 0, - scope.keyPanSpeed );\n\
				scope.update();\n\
				break;\n\
				\n\
			case scope.keys.BOTTOM:\n\
				scope.pan( 0, scope.keyPanSpeed );\n\
				scope.update();\n\
				break;\n\
				\n\
			case scope.keys.LEFT:\n\
				scope.pan( - scope.keyPanSpeed, 0 );\n\
				scope.update();\n\
				break;\n\
				\n\
			case scope.keys.RIGHT:\n\
				scope.pan( scope.keyPanSpeed, 0 );\n\
				scope.update();\n\
				break;\n\
				\n\
		}\n\
		\n\
	}\n\
	\n\
	function touchstart( event ) {\n\
	\n\
		if ( scope.enabled === false ) return;\n\
		\n\
		switch ( event.touches.length ) {\n\
		\n\
			case 1:	// one-fingered touch: rotate\n\
			\n\
				if ( scope.noRotate === true ) return;\n\
				\n\
				state = STATE.TOUCH_ROTATE;\n\
				\n\
				rotateStart.set( event.touches[ 0 ].pageX, event.touches[ 0 ].pageY );\n\
				break;\n\
				\n\
			case 2:	// two-fingered touch: dolly\n\
			\n\
				if ( scope.noZoom === true ) return;\n\
				\n\
				state = STATE.TOUCH_DOLLY;\n\
				\n\
				var dx = event.touches[ 0 ].pageX - event.touches[ 1 ].pageX;\n\
				var dy = event.touches[ 0 ].pageY - event.touches[ 1 ].pageY;\n\
				var distance = Math.sqrt( dx * dx + dy * dy );\n\
				dollyStart.set( 0, distance );\n\
				break;\n\
				\n\
			case 3: // three-fingered touch: pan\n\
			\n\
				if ( scope.noPan === true ) return;\n\
				\n\
				state = STATE.TOUCH_PAN;\n\
				\n\
				panStart.set( event.touches[ 0 ].pageX, event.touches[ 0 ].pageY );\n\
				break;\n\
				\n\
			default:\n\
			\n\
				state = STATE.NONE;\n\
				\n\
		}\n\
		\n\
		scope.dispatchEvent( startEvent );\n\
		\n\
	}\n\
	\n\
	function touchmove( event ) {\n\
	\n\
		if ( scope.enabled === false ) return;\n\
		\n\
		event.preventDefault();\n\
		event.stopPropagation();\n\
		\n\
		var element = scope.domElement === document ? scope.domElement.body : scope.domElement;\n\
		\n\
		switch ( event.touches.length ) {\n\
		\n\
			case 1: // one-fingered touch: rotate\n\
			\n\
				if ( scope.noRotate === true ) return;\n\
				if ( state !== STATE.TOUCH_ROTATE ) return;\n\
				\n\
				rotateEnd.set( event.touches[ 0 ].pageX, event.touches[ 0 ].pageY );\n\
				rotateDelta.subVectors( rotateEnd, rotateStart );\n\
				\n\
				// rotating across whole screen goes 360 degrees around\n\
				scope.rotateLeft( 2 * Math.PI * rotateDelta.x / element.clientWidth * scope.rotateSpeed );\n\
				// rotating up and down along whole screen attempts to go 360, but limited to 180\n\
				scope.rotateUp( 2 * Math.PI * rotateDelta.y / element.clientHeight * scope.rotateSpeed );\n\
				\n\
				rotateStart.copy( rotateEnd );\n\
				\n\
				scope.update();\n\
				break;\n\
				\n\
			case 2: // two-fingered touch: dolly\n\
			\n\
				if ( scope.noZoom === true ) return;\n\
				if ( state !== STATE.TOUCH_DOLLY ) return;\n\
				\n\
				var dx = event.touches[ 0 ].pageX - event.touches[ 1 ].pageX;\n\
				var dy = event.touches[ 0 ].pageY - event.touches[ 1 ].pageY;\n\
				var distance = Math.sqrt( dx * dx + dy * dy );\n\
				\n\
				dollyEnd.set( 0, distance );\n\
				dollyDelta.subVectors( dollyEnd, dollyStart );\n\
				\n\
				if ( dollyDelta.y > 0 ) {\n\
				\n\
					scope.dollyOut();\n\
					\n\
				} else {\n\
				\n\
					scope.dollyIn();\n\
					\n\
				}\n\
				\n\
				dollyStart.copy( dollyEnd );\n\
				\n\
				scope.update();\n\
				break;\n\
				\n\
			case 3: // three-fingered touch: pan\n\
			\n\
				if ( scope.noPan === true ) return;\n\
				if ( state !== STATE.TOUCH_PAN ) return;\n\
				\n\
				panEnd.set( event.touches[ 0 ].pageX, event.touches[ 0 ].pageY );\n\
				panDelta.subVectors( panEnd, panStart );\n\
				\n\
				scope.pan( panDelta.x, panDelta.y );\n\
				\n\
				panStart.copy( panEnd );\n\
				\n\
				scope.update();\n\
				break;\n\
				\n\
			default:\n\
			\n\
				state = STATE.NONE;\n\
				\n\
		}\n\
		\n\
	}\n\
	\n\
	function touchend( ) {\n\
	\n\
		if ( scope.enabled === false ) return;\n\
		\n\
		scope.dispatchEvent( endEvent );\n\
		state = STATE.NONE;\n\
		\n\
	}\n\
	\n\
	this.domElement.addEventListener( 'contextmenu', function ( event ) { event.preventDefault(); }, false );\n\
	this.domElement.addEventListener( 'mousedown', onMouseDown, false );\n\
	this.domElement.addEventListener( 'mousewheel', onMouseWheel, false );\n\
	this.domElement.addEventListener( 'DOMMouseScroll', onMouseWheel, false ); // firefox\n\
	\n\
	this.domElement.addEventListener( 'touchstart', touchstart, false );\n\
	this.domElement.addEventListener( 'touchend', touchend, false );\n\
	this.domElement.addEventListener( 'touchmove', touchmove, false );\n\
	\n\
	window.addEventListener( 'keydown', onKeyDown, false );\n\
	\n\
	// force an update at start\n\
	this.update();\n\
	\n\
	};\n\
	\n\
	THREE.OrbitControls.prototype = Object.create( THREE.EventDispatcher.prototype );",
	NULL };
