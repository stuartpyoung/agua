
dojo.provide("plugins.admin.Admin");

// DISPLAY DIFFERENT PAGES TO ALLOW THE admin AND ORDINARY
// USERS TO MODIFY THEIR SETTINGS

// DnD
dojo.require("dojo.dnd.Source"); // Source & Target
dojo.require("dojo.dnd.Moveable");
dojo.require("dojo.dnd.Mover");
dojo.require("dojo.dnd.move");

// comboBox data store
dojo.require("dojo.data.ItemFileReadStore");
dojo.require("dijit.form.ComboBox");
dojo.require("dijit.layout.ContentPane");

// rightPane buttons
dojo.require("dijit.form.Button");

// TEMPLATE MODULES
//dojo.require("plugins.admin.Access");
//dojo.require("plugins.admin.Apps");
dojo.require("plugins.admin.Settings");
//dojo.require("plugins.admin.Groups");
//dojo.require("plugins.admin.Parameter");
//dojo.require("plugins.admin.Projects");
//dojo.require("plugins.admin.GroupProjects");
//dojo.require("plugins.admin.GroupSources");
//dojo.require("plugins.admin.GroupUsers");
//dojo.require("plugins.admin.Sources");
//dojo.require("plugins.admin.Users");
dojo.require("plugins.admin.Clusters");


dojo.declare( "plugins.admin.Admin", 
	[ dijit._Widget, dijit._Templated ],
{
// PANE ID 
paneId : null,
// PANE IDENTIFIER
//paneId : "adminPane",

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "admin/templates/admin.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// PANE WIDGETS
paneWidgets : null,

// CORE WORKFLOW OBJECTS
core : new Object,

/////}}


// CONSTRUCTOR	
constructor : function(args) {

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function() {

	this.startup();
},


startup : function () {

	this.inherited(arguments);

    // ADD THIS WIDGET TO Agua.widgets
    Agua.addWidget("admin", this);

	// ADD ADMIN TAB TO TAB CONTAINER		
	Agua.tabs.addChild(this.adminTab);
	Agua.tabs.selectChild(this.adminTab);

	// RESIZE ADMIN TAB ON window.resize
	dojo.connect(window,"resize",dojo.hitch(this, function(){
		this.adminTab.resize();
	}));

	// CREATE HASH TO HOLD INSTANTIATED PANE WIDGETS
	this.paneWidgets = new Object;

	// LOAD HEADINGS FOR THIS USER
	this.headings = Agua.getHeadings();

	// LOAD LEFT PANE
	this.loadLeftPane();

	// LOAD MIDDLE PANE
	this.loadMiddlePane();

	// LOAD RIGHT PANE
	this.loadRightPane();

	for ( var tabPaneName in this.paneWidgets )
	{
	}

},


// RELOAD A WIDGET, WIDGETS IN A PANE OR ALL WIDGETS
reload : function (target) {

	if ( target == "all" )
	{
		for ( var mainPane in this.headings )
		{
			for ( var i in this.headings[mainPane] )
			{
				this.reloadWidget(this.headings[mainPane][i]);
			}
		}
	}
	else if ( target == "leftPane"
			|| target == "middlePane"
			|| target == "rightPane" )
	{
		for ( var i in this.headings[target] )
		{
			this.reloadWidget(this.headings[target][i]);
		}
	}

	// OTHERWISE, THE target MUST BE A PANE NAME
	else
	{
		try {
			this.reloadWidget(target);
		}
		catch (e) {}
	}		
},

// REINSTANTIATE A PANE WIDGET
reloadWidget : function (paneName) {

	delete this.paneWidgets[paneName];

	var adminObject = this;
	this.paneWidgets[paneName] = new plugins.admin[paneName](
		{
			parentWidget: adminObject,
			tabContainer : adminObject.leftTabContainer
		}
	);
},

// LOAD THE LEFT PANE
loadLeftPane: function () {

	if ( this.headings == null || this.headings.leftPane == null )
	{
		return;
	}

	for ( var i = 0; i < this.headings.leftPane.length; i++ )
	{
		var tabPaneName = this.headings.leftPane[i];

		var adminObject = this;
		var tabPane = new plugins["admin"][tabPaneName](
			{
				parentWidget: adminObject,
				tabContainer : adminObject.leftTabContainer
			}
		);

		// REGISTER THE NEW TAB PANE IN this.paneWidgets 
		this.paneWidgets["plugins.admin." + tabPaneName] = tabPane;
	}
}, // 	loadLeftPane 


// LOAD THE MIDDLE PANE
loadMiddlePane : function (event) {

	for ( var i = 0; i < this.headings.middlePane.length; i++ )
	{
		var tabPaneName = this.headings.middlePane[i];

		var adminObject = this;
		var tabPane = new plugins["admin"][tabPaneName](
			{
				parentWidget: adminObject,
				tabContainer : adminObject.middleTabContainer
			}
		);

		//var adminObject = this;
		//var tabPane = new plugins["admin"][tabPaneName](
		//	{
		//		parentWidget: adminObject,
		//		tabContainer : adminObject.middleTabContainer
		//	}
		//);

		// REGISTER THE NEW TAB PANE IN this.paneWidgets 
		this.paneWidgets["plugins.admin." + tabPaneName] = tabPane;
	}

},


// LOAD THE RIGHT PANE
loadRightPane : function (event) {

	if ( this.headings == null || this.headings.rightPane == null )
	{
		return;
	}

	for ( var i = 0; i < this.headings.rightPane.length; i++ )
	//for ( var i = 0; i < 3; i++ )
	{
		var tabPaneName = this.headings.rightPane[i];

		var adminObject = this;
		var tabPane = new plugins["admin"][tabPaneName](
			{
				parentWidget: adminObject,
				tabContainer : adminObject.rightTabContainer
			}
		);

		// REGISTER THE NEW TAB PANE IN this.paneWidgets 
		this.paneWidgets["plugins.admin." + tabPaneName] = tabPane;
	}
},


// LOAD CCS FILES
loadCSS : function () {

	// LOAD CSS
	var cssFiles = [ "plugins/admin/css/admin.css" ];
	for ( var i in cssFiles )
	{

		var cssFile = cssFiles[i];
		var cssNode = document.createElement('link');
		cssNode.type = 'text/css';
		cssNode.rel = 'stylesheet';
		cssNode.href = cssFile;
		document.getElementsByTagName("head")[0].appendChild(cssNode);
	}
}




}); // end of plugins.admin.Admin

