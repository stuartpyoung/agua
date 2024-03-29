dojo.provide("ImageViewAssistant");
dojo.require("dojox.mobile.app.SceneAssistant");

dojo.declare("ImageViewAssistant", dojox.mobile.app.SceneAssistant, {

  setup: function(){

	var images = [
		"images/chris1_lg.jpg",
		"images/chris2_lg.jpg",
		"images/imageHoriz.jpg",
		"images/imageVert.jpg",
		"images/square.jpg"
	];

    // Instantiate widgets in the template HTML.
    this.controller.parse();

	var viewer = this.viewer = dijit.byId("imageView");

	this.viewer.attr("leftUrl", images[0]);
	this.viewer.attr("centerUrl", images[1]);
	this.viewer.attr("rightUrl", images[2]);

	dojo.connect(dijit.byId("decZoom"), "onClick", function(){
		viewer.set("zoom", viewer.get("zoom") - 0.1);
	});
	dojo.connect(dijit.byId("incZoom"), "onClick", function(){
		viewer.set("zoom", viewer.get("zoom") + 0.1);
	});
	dojo.connect(dijit.byId("resetZoom"), "onClick", function(){
		viewer.set("zoom", 1);
	});
	dojo.connect(dijit.byId("toggleZoom"), "onClick", function(){
		if(viewer.get("zoom") > 1){
			viewer.set("animatedZoom", 1);
		}else{
			viewer.set("animatedZoom", 3);
		}
	});
	dojo.connect(dijit.byId("panLeft"), "onClick", function(){
		viewer.set("zoomCenter", {
			x: viewer.get("zoomCenterX") - 20,
			y: viewer.get("zoomCenterY")
		});
	});
	dojo.connect(dijit.byId("panRight"), "onClick", function(){
		viewer.set("zoomCenter", {
			x: viewer.get("zoomCenterX") + 20,
			y: viewer.get("zoomCenterY")
		});
	});
	dojo.connect(dijit.byId("panUp"), "onClick", function(){
		viewer.set("zoomCenter", {
			x: viewer.get("zoomCenterX"),
			y: viewer.get("zoomCenterY") - 20
		});
	});
	dojo.connect(dijit.byId("panDown"), "onClick", function(){
		viewer.set("zoomCenter", {
			x: viewer.get("zoomCenterX"),
			y: viewer.get("zoomCenterY") + 20
		});
	});

	var index = 1;

	dojo.connect(viewer, "onChange", function(direction){
		index += direction;


		if(index > 0){
			viewer.attr("leftUrl", images[index - 1]);
		}
		if(index < images.length - 1){
			viewer.attr("rightUrl", images[index + 1]);
		}
	});

  },

  activate: function(){


  }

});