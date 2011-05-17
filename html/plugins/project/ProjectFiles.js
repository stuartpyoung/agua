dojo.provide("plugins.project.ProjectFiles");

// DISPLAY THE USER'S OWN PROJECTS DIRECTORY AND ALLOW
// THE USER TO BROWSE AND MANIPULATE WORKFLOW FOLDERS
// AND FILES


// LATER FOR MENU: HOW TO DYNAMICALLY
// ENABLE / DISABLE MENU ITEM
// attr('disabled', bool) 

dojo.require("plugins.files.FileDrag");
dojo.require("dojox.data.FileStore");

// DnD
dojo.require("dojo.dnd.Source"); // Source & Target
dojo.require("dojo.dnd.Moveable");
dojo.require("dojo.dnd.Mover");
dojo.require("dojo.dnd.move");

// TOOLTIP
dojo.require("dijit.Tooltip");

// INHERITS
dojo.require("plugins.core.Common");

//HAS A
dojo.require("plugins.core.BorderContainer");
dojo.require("plugins.core.ExpandoPane");

// INHERITS
dojo.require("plugins.core.Common");

// MENUS
dojo.require("plugins.files.FileMenu");
dojo.require("plugins.files.FolderMenu");
dojo.require("plugins.files.WorkflowMenu");

// HAS A TITLE PANE IN ITS TEMPLATE.
// ALSO INSERTS TITLE PANES INTO this.rowsNode
dojo.require("plugins.files.TitlePane");


dojo.declare( "plugins.project.ProjectFiles",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//////}}

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "project/templates/filesystem.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// PARENT WIDGET
parentWidget : null,

// PROJECT NAME 
project : null,

// DEFAULT TIME (milliseconds) TO SLEEP BETWEEN FILESYSTEM LOADS
sleep : 300,

// CSS FILES
cssFiles : [
	"plugins/project/css/dialog.css",
	"dojo-1.5.0/dojox/widget/Dialog/Dialog.css",
	"plugins/files/FileDrag/FileDrag.css"
],

// CONSTRUCTOR	
constructor : function(args) {

	// SET PANE ID
	this.paneId = args.paneId;

	// SET ARGS
	this.open = args.open;
	if ( this.open == null )	this.open = false;
	this.parentWidget = args.parentWidget;
	this.attachNode = args.attachNode;

	// LOAD CSS
	this.loadCSS();
},


// DO STARTUP
postCreate: function() {

	// SET RELOAD CALLBACK
	this.titlePane.reloadCallback = dojo.hitch(this, "reload");

	this.startup();
},

// START MENUS
startup : function () {

	// SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
	this.inherited(arguments);

	// ADD THE PANE TO THE TAB CONTAINER
	this.attachNode.appendChild(this.mainNode);

	// LOAD APPLICATIONS DND SOURCE AND TARGET INTO FRAMEWORK
	this.load();
},


// RELOAD THE WHOLE PROJECT TAB
reload : function () {
	// REMOVE EXISTING DIRECTORIES
	while ( this.rowsNode.firstChild )
	{
		this.rowsNode.removeChild(this.rowsNode.firstChild);
	}

	// LOAD DIRECTORIES ANEW
	this.load();
},


// GET DIRECTORIES TO SEARCH FOR FILES
getDirectories : function () {

	var projects = Agua.getProjects();

	return projects;
},


// LOAD FILESYSTEMS
load : function () {

	var directories = this.getDirectories();
	if ( directories == null )
	{
		return;
	}


	var sleep = this.sleep;
	var sequenceArray = new Array;
	for ( var i = 0; i < directories.length; i++ )
	{
		var directory = dojo.clone(directories[i]);
		sequenceArray.push({ func: [dojo.hitch(this,"setFileDrag", directory)], pauseAfter: 150 });
		sequenceArray.push({ func: [showMessage, window, "ProjectFiles.setFileDrag    Finished " + sleep + " milliseconds pause after load " + directory.name], pauseAfter: sleep });
	}

	function showMessage(msg) {
	}

	/* SET TIMED SEQUENCE: runWorkflow NEEDS A DELAY TO COMPLETE
	    EXECUTION (ON WINDOWS) BEFORE runStatus IS CALLED	
		var seq = [
			{ func: [showMessage, window, "repeat 10 times and pause 100ms after"], repeat: 2, pauseAfter: 100 },
			{ func: dojo.hitch(workflow, "runWorkflow"), pauseAfter: 8000 },
			{ func: [showMessage, window, "after 1000ms pause this should be seen"], pauseAfter: 1000 },
			{ func: dojo.hitch(workflow, "runStatus"), pauseAfter: 10000, repeat:100 },
			{ func: returnWhenDone }  no array, just a function to call 
			{ func: [dojo.hitch(console,"log","hoot"), window, "hoot2?"], pauseAfter: 1000, repeat:20 }
		];


		REPORTER FUNCTION

		function returnWhenDone() {
			//logMsg("in returnWhenDone");
			//setTimeout(continueSequence, 1000);
			return false;
		}


	*/

	var sequenceObject = null,
	sequenceObject = new dojox.timing.Sequence({});

},


// SET THE FILE SYSTEM PANE
setFileDrag : function (directory) {

	var owner = directory.username;
	var name = directory.name;
	var description = directory.description;
	if ( ! description ) { description = '' };

	var titlePane = this.createTitlePane(
	{
		owner: owner,
		name: name,
		description: description,
		open: this.open,
		directory: directory
	});

	var fileDrag = this.createFileDrag(directory);

	// ADD fileDrag TO TITLE PANE
	titlePane.containerNode.appendChild(fileDrag.domNode);	
},




// SET THE FILE SYSTEM PANE
createFileDrag : function (directory) {

	var owner = directory.username;
	var name = directory.name;
	var description = directory.description;
	if ( ! description ) { description = '' };

	var thisStore = this.createStore(directory);

	// GENERATE NEW FileDrag OBJECT
	var fileDrag = new plugins.files.FileDrag(
		{
			style: "height: auto; width: 100%; minHeight: 50px;",
			store: thisStore,
			fileMenu: this.fileMenu,
			folderMenu: this.folderMenu,
			workflowMenu: this.workflowMenu,
			owner: owner,
			parentWidget: this
		}
	);

	// SET PATH FOR THIS PROJECT
	fileDrag.path = name;                    

	// START UP FileDrag
	fileDrag.startup();

	return fileDrag;
},




// CREATE TITLE PANE
createTitlePane : function (directory) {

	var titlePaneNode = document.createElement('div');
	this.rowsNode.appendChild(titlePaneNode);

	// CREATE TITLE PANE			
	var titlePane = new plugins.files.TitlePane(
		{
			owner 	: directory.owner,
			type 	: directory.type,
			name	: directory.name,
			description: directory.description,

			open: this.open,
			//reloadCallback : callback,
			directory : directory,
			parentWidget : this
		},
		titlePaneNode
	);

	return titlePane;
},


// CREATE URL FOR STORE
// **** OVERRIDE THIS IN SUBCLASS ****
createUrl : function (directory) {
	var url = Agua.cgiUrl;
	url += "project.cgi?";
	url += "mode=fileSystem";
	url += "&username=" + Agua.cookie('username');
	url += "&sessionId=" + Agua.cookie('sessionId');

	return url;
},


// CREATE STORE FOR FILE DRAG
createStore : function (directory) {

	// SET URL
	var url = this.createUrl(directory);

	// CREATE STORE
	var thisStore = new dojox.data.FileStore(
		{
			//id: paneNodeId + "-fileStore",
			url: url,
			pathAsQueryParam: true
		}
	);

	// SET FILE STORE path TO project
	thisStore.preamble = function()
	{
		this.store.path = this.arguments[0].path;                        
	};

	return thisStore;		
}



});
