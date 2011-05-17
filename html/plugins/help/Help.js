
dojo.provide("plugins.help.Help");

// DISPLAY HELP PAGES


// Buttons
dojo.require("dijit.form.Button");

// INHERITS
dojo.require("plugins.core.Common");

dojo.declare( "plugins.help.Help", 
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	// PANE ID 
	paneId : null,
	// PANE IDENTIFIER
	//paneId : "helpPane",

	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "help/templates/help.html"),

	cssFiles : [ "plugins/help/css/help.css"],

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// PANE WIDGETS
	paneWidgets : null,

 	// CONSTRUCTOR	
	constructor : function(args) {

        // LOAD SORIA AND FILEPICKER CSS
        this.loadCSS();		
	},

	postCreate : function()
	{

		this.startup();
	},


	startup : function ()
	{

		this.inherited(arguments);

		// ADD ADMIN TAB TO TAB CONTAINER		
		Agua.tabs.addChild(this.mainTab);
		Agua.tabs.selectChild(this.mainTab);

		// RESIZE ADMIN TAB ON window.resize
		dojo.connect(window,"resize",dojo.hitch(this, function(){
			this.mainTab.resize();
		}));

		// LOAD PANE 
		//this.loadPane();
	},

	faq : function(event)
	{

		window.open('http://gouda.ccs.miami.edu:8090/display/Agua/FAQ');
	},

	gettingStarted : function (event)
	{

		window.open('http://gouda.ccs.miami.edu:8090/display/Agua/Creating+workflows');
	},

	sharing : function (event)
	{

		window.open('http://gouda.ccs.miami.edu:8090/display/Agua/Managing+projects');
	},

	openWorkflows : function (event)
	{

		if ( Agua.workflowController )
		{
			Agua.workflowController.createTab();
		}
		else
		{
			var uniqueId = dijit.getUniqueId("workflowTab");
			Agua.workflowController = new plugins.workflow.Controller( { id: uniqueId, target: 'tabs', deleteRoot: false } );
			Agua.workflowController.createTab();
		}
	},

	openProjects : function (event)
	{

		if ( Agua.projectController )
		{
			Agua.projectController.createTab();
		}
		else
		{
			var uniqueId = dijit.getUniqueId("projectTab");
			Agua.projectController = new plugins.project.Controller( { id: uniqueId, target: 'tabs', deleteRoot: false } );
			Agua.projectController.createTab({ id: uniqueId, target: 'tabs', deleteRoot: false });
		}	
	}

}); // end of plugins.help.Help
