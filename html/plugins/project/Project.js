
dojo.provide("plugins.project.Project");

// DISPLAY THE USER'S PROJECTS DIRECTORY AND ALLOW
// THE USER TO BROWSE FILES AND MANIPULATE WORKFLOW
// FOLDERS AND FILES
dojo.require("plugins.files.FileDrag");
dojo.require("plugins.files._FileInfoPane");
dojo.require("dojox.data.FileStore");

// DnD
dojo.require("dojo.dnd.Source"); // Source & Target
dojo.require("dojo.dnd.Moveable");
dojo.require("dojo.dnd.Mover");
dojo.require("dojo.dnd.move");

// TITLE PANE
dojo.require("dijit.TitlePane");

// TOOLTIP
dojo.require("dijit.Tooltip");

// FILE UPLOAD
//dojo.require("plugins.project.FileInput");	
//dojo.require("plugins.project.FileInputAuto");	
dojo.require("dojox.layout.ResizeHandle");

// INHERITS
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.project.ProjectFiles");
dojo.require("plugins.project.SharedProjectFiles");
dojo.require("plugins.project.SourceFiles");
dojo.require("plugins.project.SharedSourceFiles");

//HAS A
//dojo.require("plugins.core.BorderContainer");
dojo.require("dijit.layout.ContentPane");
dojo.require("plugins.core.ExpandoPane");
dojo.require("plugins.core.Common");

// MENUS
dojo.require("plugins.files.FileMenu");
dojo.require("plugins.files.FolderMenu");
dojo.require("plugins.files.WorkflowMenu");


dojo.declare( "plugins.project.Project",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//////}}

// PANE IDENTIFIER
paneId : '',

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "project/templates/project.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// PARENT NODE, I.E., TABS NODE
parentWidget : null,

// PROJECT NAME AND WORKFLOW NAME IF AVAILABLE
project : null,
project : null,

// POLL SERVER FOR WORKFLOW STATUS
polling : false,

// INSERT TEXT BREAKS WIDTH, CORRESPONDS TO CSS WIDTH OF INPUT 'value' TABLE ELEMENT
textBreakWidth : 22,

// CSS FILES
cssFiles : [
	"plugins/project/css/project.css"
],

constructor : function(args) {		

	// SET PANE ID
	this.paneId = args.paneId;

	// SET ARGS
	this.parentWidget = Agua.tabs;
	this.project = args.project;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();

},

postMixInProperties: function() {
},

postCreate: function() {
	this.startup();
},

startup : function () {

	// SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
	this.inherited(arguments);

	// ADD THE PANE TO THE TAB CONTAINER
	this.parentWidget.addChild(this.mainTab);
	this.parentWidget.selectChild(this.mainTab);

	// SET MENUS
	this.setMenus();

	// LOAD APPLICATIONS DND SOURCE AND TARGET INTO FRAMEWORK
	this.loadProjectTab();
},


setMenus : function () {

	this.fileMenu = new plugins.files.FileMenu();
	this.folderMenu = new plugins.files.FolderMenu();
	this.workflowMenu = new plugins.files.WorkflowMenu();


	// STOP PROPAGATION TO NORMAL RIGHTCLICK CONTEXT MENU
	dojo.connect(this.fileMenu.menu.domNode, "oncontextmenu", function (event)
	{
		event.stopPropagation();
	});

	//// STOP PROPAGATION TO NORMAL RIGHTCLICK CONTEXT MENU
	//dojo.connect(this.folderMenu.menu.domNode, "oncontextmenu", function (event)
	//{
	//	event.stopPropagation();
	//});
	//
	//// STOP PROPAGATION TO NORMAL RIGHTCLICK CONTEXT MENU
	//dojo.connect(this.workflowMenu.menu.domNode, "oncontextmenu", function (event)
	//{
	//	event.stopPropagation();
	//});
},

// RELOAD THE WHOLE PROJECT TAB
reload : function () {
	this.loadProjectTab();
},

// LOAD THE PROJECT TAB
loadProjectTab: function ()	{

	// REMOVE EXISTING TITLE PANES
	while ( this.projectsNode.firstChild )
	{
		this.projectsNode.removeChild(this.projectsNode.firstChild);
	}

	// LOAD FILE SYSTEMS
	this.loadProjects();


	//this.loadSharedProjects();
	//this.loadSources();
	//this.loadSharedSources();



},



// LOAD PROJECTS
loadProjects : function () {
	this.projectFiles = new plugins.project.ProjectFiles(
	{
		open: true,
		title: 'Projects',
		type: 'Project',
		attachNode : this.projectsNode,
		parentWidget : this,
		fileMenu: this.fileMenu,
		folderMenu: this.folderMenu,
		workflowMenu: this.workflowMenu
	});

	return;
},



loadSharedProjects : function () {
	this.sharedProjectFiles = new plugins.project.SharedProjectFiles(
	{
		open: false,
		title: 'Shared Projects',
		type: 'Shared Project',
		attachNode : this.projectsNode,
		parentWidget : this,
		fileMenu: this.fileMenu,
		folderMenu: this.folderMenu,
		workflowMenu: this.workflowMenu
	});

	return;
},

// LOAD SOURCE FILE PANES
loadSources : function () {
	this.sourceFiles = new plugins.project.SourceFiles(
	{
		open: false,
		title: 'Sources',
		type: 'Source',
		attachNode : this.projectsNode,
		parentWidget : this,
		project : this.project,
		fileMenu: this.fileMenu,
		folderMenu: this.folderMenu,
		workflowMenu: this.workflowMenu
	});

	return;
},

loadSharedSources : function () {
	this.sharedSourceFiles = new plugins.project.SharedSourceFiles(
	{
		open: false,
		title: 'SharedSources',
		type: 'Shared Source',
		attachNode : this.projectsNode,
		parentWidget : this,
		project : this.project,
		fileMenu: this.fileMenu,
		folderMenu: this.folderMenu,
		workflowMenu: this.workflowMenu
	});

}





/*
//
//		// ***********************
//		// *** SHARED SOURCES ****
//		// ***********************
//		var sharedSources = Agua.getSharedSources();
//		if ( sharedSources == null ) 	return;
//		
//		
//		for ( var i = 0; i < sharedSources.length; i++ )
//		{
//			var source = sharedSources[i];
//
//			// CREATE TITLE 
//			var titlePaneNode = document.createElement('div');
//			this.projectsNode.appendChild(titlePaneNode);
//
//			var owner = sharedSources[i].owner;
//			var groupname = sharedSources[i].groupname;
//			var name = sharedSources[i].name;
//			var description = sharedSources[i].description;
//			
//			// CREATE TITLE TABLE
//			var title = "<table class='titlePane'><tr><td class='owner'>" + owner + "</td><td class='name'> Shared Source: " + name + "</td><td class='description'>" + description + "</td></tr></table>";
//
//			// CREATE TITLE PANE			
//			var titlePane = new dijit.TitlePane(
//				{
//					title: title,
//					//duration : 100, 
//					open: false
//				},
//				titlePaneNode
//			);
//			
//			// SET URL
//			var url = Agua.cgiUrl;
//			url += "project.cgi?";
//			
//			// SET QUERY
//			var query = "mode=fileSystem";
//			query += "&username=" + owner;
//			query += "&requestor=" + Agua.cookie('username');
//			query += "&sessionId=" + Agua.cookie('sessionId');
//
//			// SET LOCATION TO BE USED INSTEAD OF FILEROOT
//			// ASSUMES ALL SOURCES ARE AT LEAST TWO FOLDERS DEEP IN THE FILESYSTEM
//			var sourceBase = source.location.match(/^(.+)\/[^\/]+$/)[1];
//			var sourceFolder = source.location.match(/^.+\/([^\/]+)$/)[1];
//			query += "&location=" + sourceBase;
//			
//			// GROUPNAME
//			query += "&groupname=" + groupname;
//
//			
//			// CREATE STORE
//			var projectStore = new dojox.data.FileStore(
//				{
//					url: url + query,
//					pathAsQueryParam: true
//				}
//			);
//
//			// GENERATE NEW FileDrag OBJECT
//			var sourceDragId = dojo.dnd.getUniqueId();
//			var sourceDrag = new plugins.files.FileDrag(
//				{
//					id: sourceDragId,
//					style: "height: auto; width: 100%; minHeight: 50px;",
//					store: projectStore,
//					fileMenu: this.fileMenu,
//					folderMenu: this.folderMenu,
//					workflowMenu: this.workflowMenu,
//					owner: owner,
//					groupname: groupname,
//					parentWidget: this
//				}
//			);
//
//
//
//
//
//
//			// SET PATH FOR THIS PROJECT
//			sourceDrag.path = sourceFolder;                    
//
//return;
//		
//			// START UP FileDrag
//			sourceDrag.startup();
//
//
//		//console.dir(projectDrag);
//
//							
//			// ADD sourceDrag TO TITLE PANE
//			titlePane.containerNode.appendChild(sourceDrag.domNode);
//
//
//		} // for loop on sharedSources.length		
//
//
//	}
//
*/

}); // end of plugins.project.Project
