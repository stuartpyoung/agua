/*
	Copyright (c) 2004-2010, The Dojo Foundation All Rights Reserved.
	Available via Academic Free License >= 2.1 OR the modified BSD license.
	see: http://dojotoolkit.org/license for details
*/


dojo.provide("plugins.core.Controls");

// GENERATE THE FOLLOWING ELEMENTS AT THE TOP OF THE PAGE:
//
// 	1. PLUGINS TOOLBAR - TO BE POPULATED WITH PLUGIN BUTTONS
//
// 	2. TABS CONTAINER - TO BE POPULATED BY PLUGIN TAB PANES
// 
//	3. LOGIN STATUS - TO BE POPULATED BY LOGIN PLUGIN
//

// REQUIRE WIDGETS
dojo.require("dojox.layout.ResizeHandle");
dojo.require("dijit.layout.LayoutContainer");
dojo.require("dijit.Toolbar");
dojo.require("dijit.layout.TabContainer");

dojo.declare( "plugins.core.Controls",
	null,
{
	// PANE IDENTIFIER
	paneId : 1,

	//Path to the template of this widget. 
	//templatePath: dojo.moduleUrl("plugins", "core/templates/controls.html"),

	// Calls dijit._Templated.widgetsInTemplate
	//widgetsInTemplate : true,

	//You can override this method to manipulate widget once it is
	//placed in the UI, but be warned that any child widgets contained
	//in it may not be ready yet.        
	postCreate: function()
	{

		this.startup();
	},


	// CONSTRUCTOR	
	constructor : function(args)
	{		

		// SET PANE ID
		this.paneId = args.paneId;
	},


	// DO inherited, LOAD CSS AND
	// SET Agua.controls, Agua.toolbar AND Agua.tabs 
	startup : function (args)
	{

return;


		// SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
		this.inherited(arguments);

		this.loadCSS();

		// SET CONTROLS
		//Agua.controls = this.controls;

		// SET TOOLBAR
		//Agua.toolbar = this.toolbar;

		// SET TABS
		//Agua.tabs = dijit.byId("mainTabContainer");

		// FIX ERROR: TABSLIDER TAKES TWO LINES
		var tabSlider = dijit.byId('dijit_layout_TabContainer_0_tablist_menuBtn');
		tabSlider.domNode.setAttribute('style', 'display: none');

		// ATTACH TO THIS TEMPLATE TO controls DIV ON HTML PAGE
		//Agua.controls = dojo.byId("controls");
		Agua.controls.appendChild(this.containerNode);

		// LOAD FIRST BUTTON IN TOOLBAR
		var label = "Agua";
		this.aguaButton = Agua.addToolbarButton(label);

		// SET BUTTON LISTENER
		var listener = dojo.connect(this.aguaButton, "onClick", this, "reload");
	},


	// RELOAD AGUA
	reload : function ()
	{

		var url = window.location;
		window.open(location, '_blank', 'toolbar=1,location=0,directories=0,status=0,menubar=1,scrollbars=1,resizable=1,navigation=0'); 

		//window.location.reload();
	},


	loadCSS : function ()
    {

		// LOAD CSS
		var cssFiles = [ "plugins/core/css/controls.css" ];
		for ( var i in cssFiles )
		{
			var cssFile = cssFiles[i];
			var cssNode = document.createElement('link');
			cssNode.type = 'text/css';
			cssNode.rel = 'stylesheet';
			cssNode.href = cssFile;
			cssNode.media = 'screen';
			cssNode.title = 'coreCSS';
			document.getElementsByTagName("head")[0].appendChild(cssNode);
		}
    }

}); // plugins.core.Controls

dojo.addOnLoad(
	function()
	{
		var controls = new plugins.core.Controls( { id: 1, target: 'controls' } );
	}
);

