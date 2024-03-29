dojo.provide("dojox.gfx.path");

dojo.require("dojox.gfx.matrix");
dojo.require("dojox.gfx.shape");

dojo.declare("dojox.gfx.path.Path", dojox.gfx.Shape, {
	// summary: a generalized path shape

	constructor: function(rawNode){
		// summary: a path constructor
		// rawNode: Node: a DOM node to be used by this path object
		this.shape = dojo.clone(dojox.gfx.defaultPath);
		this.segments = [];
		this.tbbox = null;
		this.absolute = true;
		this.last = {};
		this.rawNode = rawNode;
		this.segmented = false;
	},

	// mode manipulations
	setAbsoluteMode: function(mode){
		// summary: sets an absolute or relative mode for path points
		// mode: Boolean: true/false or "absolute"/"relative" to specify the mode
		this._confirmSegmented();
		this.absolute = typeof mode == "string" ? (mode == "absolute") : mode;
		return this; // self
	},
	getAbsoluteMode: function(){
		// summary: returns a current value of the absolute mode
		this._confirmSegmented();
		return this.absolute; // Boolean
	},

	getBoundingBox: function(){
		// summary: returns the bounding box {x, y, width, height} or null
		this._confirmSegmented();
		return (this.bbox && ("l" in this.bbox)) ? {x: this.bbox.l, y: this.bbox.t, width: this.bbox.r - this.bbox.l, height: this.bbox.b - this.bbox.t} : null; // dojox.gfx.Rectangle
	},

	_getRealBBox: function(){
		// summary: returns an array of four points or null
		//	four points represent four corners of the untransformed bounding box
		this._confirmSegmented();
		if(this.tbbox){
			return this.tbbox;	// Array
		}
		var bbox = this.bbox, matrix = this._getRealMatrix();
		this.bbox = null;
		for(var i = 0, len = this.segments.length; i < len; ++i){
			this._updateWithSegment(this.segments[i], matrix);
		}
		var t = this.bbox;
		this.bbox = bbox;
		this.tbbox = t ? [
			{x: t.l, y: t.t},
			{x: t.r, y: t.t},
			{x: t.r, y: t.b},
			{x: t.l, y: t.b}
		] : null;
		return this.tbbox;	// Array
	},

	getLastPosition: function(){
		// summary: returns the last point in the path, or null
		this._confirmSegmented();
		return "x" in this.last ? this.last : null; // Object
	},

	_applyTransform: function(){
		this.tbbox = null;
		return dojox.gfx.Shape.prototype._applyTransform.call(this);
	},

	// segment interpretation
	_updateBBox: function(x, y, matrix){
		// summary: updates the bounding box of path with new point
		// x: Number: an x coordinate
		// y: Number: a y coordinate

		if(matrix){
			var t = dojox.gfx.matrix.multiplyPoint(matrix, x, y);
			x = t.x;
			y = t.y;
		}

		// we use {l, b, r, t} representation of a bbox
		if(this.bbox && ("l" in this.bbox)){
			if(this.bbox.l > x) this.bbox.l = x;
			if(this.bbox.r < x) this.bbox.r = x;
			if(this.bbox.t > y) this.bbox.t = y;
			if(this.bbox.b < y) this.bbox.b = y;
		}else{
			this.bbox = {l: x, b: y, r: x, t: y};
		}
	},
	_updateWithSegment: function(segment, matrix){
		// summary: updates the bounding box of path with new segment
		// segment: Object: a segment
		var n = segment.args, l = n.length;
		// update internal variables: bbox, absolute, last
		switch(segment.action){
			case "M":
			case "L":
			case "C":
			case "S":
			case "Q":
			case "T":
				for(var i = 0; i < l; i += 2){
					this._updateBBox(n[i], n[i + 1], matrix);
				}
				this.last.x = n[l - 2];
				this.last.y = n[l - 1];
				this.absolute = true;
				break;
			case "H":
				for(var i = 0; i < l; ++i){
					this._updateBBox(n[i], this.last.y, matrix);
				}
				this.last.x = n[l - 1];
				this.absolute = true;
				break;
			case "V":
				for(var i = 0; i < l; ++i){
					this._updateBBox(this.last.x, n[i], matrix);
				}
				this.last.y = n[l - 1];
				this.absolute = true;
				break;
			case "m":
				var start = 0;
				if(!("x" in this.last)){
					this._updateBBox(this.last.x = n[0], this.last.y = n[1], matrix);
					start = 2;
				}
				for(var i = start; i < l; i += 2){
					this._updateBBox(this.last.x += n[i], this.last.y += n[i + 1], matrix);
				}
				this.absolute = false;
				break;
			case "l":
			case "t":
				for(var i = 0; i < l; i += 2){
					this._updateBBox(this.last.x += n[i], this.last.y += n[i + 1], matrix);
				}
				this.absolute = false;
				break;
			case "h":
				for(var i = 0; i < l; ++i){
					this._updateBBox(this.last.x += n[i], this.last.y, matrix);
				}
				this.absolute = false;
				break;
			case "v":
				for(var i = 0; i < l; ++i){
					this._updateBBox(this.last.x, this.last.y += n[i], matrix);
				}
				this.absolute = false;
				break;
			case "c":
				for(var i = 0; i < l; i += 6){
					this._updateBBox(this.last.x + n[i], this.last.y + n[i + 1], matrix);
					this._updateBBox(this.last.x + n[i + 2], this.last.y + n[i + 3], matrix);
					this._updateBBox(this.last.x += n[i + 4], this.last.y += n[i + 5], matrix);
				}
				this.absolute = false;
				break;
			case "s":
			case "q":
				for(var i = 0; i < l; i += 4){
					this._updateBBox(this.last.x + n[i], this.last.y + n[i + 1], matrix);
					this._updateBBox(this.last.x += n[i + 2], this.last.y += n[i + 3], matrix);
				}
				this.absolute = false;
				break;
			case "A":
				for(var i = 0; i < l; i += 7){
					this._updateBBox(n[i + 5], n[i + 6], matrix);
				}
				this.last.x = n[l - 2];
				this.last.y = n[l - 1];
				this.absolute = true;
				break;
			case "a":
				for(var i = 0; i < l; i += 7){
					this._updateBBox(this.last.x += n[i + 5], this.last.y += n[i + 6], matrix);
				}
				this.absolute = false;
				break;
		}
		// add an SVG path segment
		var path = [segment.action];
		for(var i = 0; i < l; ++i){
			path.push(dojox.gfx.formatNumber(n[i], true));
		}
		if(typeof this.shape.path == "string"){
			this.shape.path += path.join("");
		}else{
			Array.prototype.push.apply(this.shape.path, path);
		}
	},

	// a dictionary, which maps segment type codes to a number of their arguments
	_validSegments: {m: 2, l: 2, h: 1, v: 1, c: 6, s: 4, q: 4, t: 2, a: 7, z: 0},

	_pushSegment: function(action, args){
		// summary: adds a segment
		// action: String: valid SVG code for a segment's type
		// args: Array: a list of parameters for this segment
		this.tbbox = null;
		var group = this._validSegments[action.toLowerCase()];
		if(typeof group == "number"){
			if(group){
				if(args.length >= group){
					var segment = {action: action, args: args.slice(0, args.length - args.length % group)};
					this.segments.push(segment);
					this._updateWithSegment(segment);
				}
			}else{
				var segment = {action: action, args: []};
				this.segments.push(segment);
				this._updateWithSegment(segment);
			}
		}
	},

	_collectArgs: function(array, args){
		// summary: converts an array of arguments to plain numeric values
		// array: Array: an output argument (array of numbers)
		// args: Array: an input argument (can be values of Boolean, Number, dojox.gfx.Point, or an embedded array of them)
		for(var i = 0; i < args.length; ++i){
			var t = args[i];
			if(typeof t == "boolean"){
				array.push(t ? 1 : 0);
			}else if(typeof t == "number"){
				array.push(t);
			}else if(t instanceof Array){
				this._collectArgs(array, t);
			}else if("x" in t && "y" in t){
				array.push(t.x, t.y);
			}
		}
	},

	// segments
	moveTo: function(){
		// summary: formes a move segment
		this._confirmSegmented();
		var args = [];
		this._collectArgs(args, arguments);
		this._pushSegment(this.absolute ? "M" : "m", args);
		return this; // self
	},
	lineTo: function(){
		// summary: formes a line segment
		this._confirmSegmented();
		var args = [];
		this._collectArgs(args, arguments);
		this._pushSegment(this.absolute ? "L" : "l", args);
		return this; // self
	},
	hLineTo: function(){
		// summary: formes a horizontal line segment
		this._confirmSegmented();
		var args = [];
		this._collectArgs(args, arguments);
		this._pushSegment(this.absolute ? "H" : "h", args);
		return this; // self
	},
	vLineTo: function(){
		// summary: formes a vertical line segment
		this._confirmSegmented();
		var args = [];
		this._collectArgs(args, arguments);
		this._pushSegment(this.absolute ? "V" : "v", args);
		return this; // self
	},
	curveTo: function(){
		// summary: formes a curve segment
		this._confirmSegmented();
		var args = [];
		this._collectArgs(args, arguments);
		this._pushSegment(this.absolute ? "C" : "c", args);
		return this; // self
	},
	smoothCurveTo: function(){
		// summary: formes a smooth curve segment
		this._confirmSegmented();
		var args = [];
		this._collectArgs(args, arguments);
		this._pushSegment(this.absolute ? "S" : "s", args);
		return this; // self
	},
	qCurveTo: function(){
		// summary: formes a quadratic curve segment
		this._confirmSegmented();
		var args = [];
		this._collectArgs(args, arguments);
		this._pushSegment(this.absolute ? "Q" : "q", args);
		return this; // self
	},
	qSmoothCurveTo: function(){
		// summary: formes a quadratic smooth curve segment
		this._confirmSegmented();
		var args = [];
		this._collectArgs(args, arguments);
		this._pushSegment(this.absolute ? "T" : "t", args);
		return this; // self
	},
	arcTo: function(){
		// summary: formes an elliptic arc segment
		this._confirmSegmented();
		var args = [];
		this._collectArgs(args, arguments);
		this._pushSegment(this.absolute ? "A" : "a", args);
		return this; // self
	},
	closePath: function(){
		// summary: closes a path
		this._confirmSegmented();
		this._pushSegment("Z", []);
		return this; // self
	},

	_confirmSegmented: function() {
		if (!this.segmented) {
			var path = this.shape.path;
			// switch to non-updating version of path building
			this.shape.path = [];
			this._setPath(path);
			// switch back to the string path
			this.shape.path = this.shape.path.join("");
			// become segmented
			this.segmented = true;
		}
	},

	// setShape
	_setPath: function(path){
		// summary: forms a path using an SVG path string
		// path: String: an SVG path string
		var p = dojo.isArray(path) ? path : path.match(dojox.gfx.pathSvgRegExp);
		this.segments = [];
		this.absolute = true;
		this.bbox = {};
		this.last = {};
		if(!p) return;
		// create segments
		var action = "",	// current action
			args = [],		// current arguments
			l = p.length;
		for(var i = 0; i < l; ++i){
			var t = p[i], x = parseFloat(t);
			if(isNaN(x)){
				if(action){
					this._pushSegment(action, args);
				}
				args = [];
				action = t;
			}else{
				args.push(x);
			}
		}
		this._pushSegment(action, args);
	},
	setShape: function(newShape){
		// summary: forms a path using a shape
		// newShape: Object: an SVG path string or a path object (see dojox.gfx.defaultPath)
		dojox.gfx.Shape.prototype.setShape.call(this, typeof newShape == "string" ? {path: newShape} : newShape);

		this.segmented = false;
		this.segments = [];
		if(!dojox.gfx.lazyPathSegmentation){
			this._confirmSegmented();
		}
		return this; // self
	},

	// useful constant for descendants
	_2PI: Math.PI * 2
});

dojo.declare("dojox.gfx.path.TextPath", dojox.gfx.path.Path, {
	// summary: a generalized TextPath shape

	constructor: function(rawNode){
		// summary: a TextPath shape constructor
		// rawNode: Node: a DOM node to be used by this TextPath object
		if(!("text" in this)){
			this.text = dojo.clone(dojox.gfx.defaultTextPath);
		}
		if(!("fontStyle" in this)){
			this.fontStyle = dojo.clone(dojox.gfx.defaultFont);
		}
	},
	getText: function(){
		// summary: returns the current text object or null
		return this.text;	// Object
	},
	setText: function(newText){
		// summary: sets a text to be drawn along the path
		this.text = dojox.gfx.makeParameters(this.text,
			typeof newText == "string" ? {text: newText} : newText);
		this._setText();
		return this;	// self
	},
	getFont: function(){
		// summary: returns the current font object or null
		return this.fontStyle;	// Object
	},
	setFont: function(newFont){
		// summary: sets a font for text
		this.fontStyle = typeof newFont == "string" ?
			dojox.gfx.splitFontString(newFont) :
			dojox.gfx.makeParameters(dojox.gfx.defaultFont, newFont);
		this._setFont();
		return this;	// self
	}
});
