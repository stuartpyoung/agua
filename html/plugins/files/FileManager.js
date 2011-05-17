dojo.provide("plugins.files.FileManager");

// DISPLAY THE USER'S PROJECTS DIRECTORY AND ALLOW
// THE USER TO BROWSE FILES AND MANIPULATE WORKFLOW
// FOLDERS AND FILES

dojo.require("plugins.files.FileSelector");
dojo.require("dojox.data.FileStore");

dojo.require("dojox.widget.Dialog");


// DnD
dojo.require("dojo.dnd.Source"); // Source & Target
dojo.require("dojo.dnd.Moveable");
dojo.require("dojo.dnd.Mover");
dojo.require("dojo.dnd.move");

// FILE UPLOAD
dojo.require("dojox.widget.FileInput");	
dojo.require("dojox.widget.FileInputAuto");	

// INHERITS
dojo.require("plugins.core.Common");

// WidgetFramework ancestor
dojo.require("plugins.core.WidgetFramework");

dojo.declare( "plugins.files.FileManager",
	[ plugins.files.Project, plugins.core.Common ],
{
	// PANE IDENTIFIER
	paneId : '',

	// STORE DND SOURCE ID AND DND TARGET ID
	sourceId: '',

	// CSS FILES
	ccsFiles : [
		"plugins/files/FileDrag/FileDrag.css",
		"dojo-1.5.0/dijit/themes/soria/soria.css",
		"plugins/project/FileSelector/FileManager.css"
	],


	//// LAYOUT TO GENERATE WORKFLOW PANE
	//layout :
	//{
	//	widgetType: "dojox.widget.Dialog",
	//	//params: {id: "tabContainer", closable: "true", title: this.paneId, style:'height:450px; width:850px;'},
	//	params: { id: "fileDialog", draggable: 'true', sizeDuration: 200, sizeMethod: 'combine', viewportPadding: '125', showTitle: 'true', title: 'File Manager', dimensions: [ 800,800 ] },
	//	innerHTML: '',
	//	children:
	//	[
	//		{
	//			widgetType: "dijit.layout.ContentPane",
	//			params: {id: "leftPane", sizeShare:  100, style: "background: #EEEEFF;" },
	//			style: "background: #EEEEFF;",
	//			innerHTML: "Project folders"
	//		}
	//	]
	//},


	// DIALOGUE TO DISPLAY FILE MANAGER
	fileManagerDialog : null,

	// callback FUNCTION AND DATA FROM OBJECT THAT GENERATED THE FileManager
	callback : null,


	preamble: function(){

		this.callback = arguments[0].callback;
	},

 	// CONSTRUCTOR	
	constructor : function(args) {

		this.callback = args.callback;

		if ( ! args.paneId )
		{
			return;
		}

		// SET PANE ID
		this.paneId = args.paneId;


		// SET IDS IN LAYOUT USING PANE ID
		//this.setIds();

		// SET TAB CONTAINER NODE ID ( THE SAME AS THE ID SET ABOVE IN this.setIds() )
		//var tabContainerNodeId = this.paneId + "tabContainer";

		// CREATE FRAMEWORK USING INHERITED create FROM WidgetFramework
		//this.create(this.layout, { root: tabContainerNodeId, target: 'tab' + this.paneId } );

		// APPEND TAB TO TABS CONTAINER
		//var tabsNodeId = args.tabsNodeId;
		//var tabsNode = dijit.byId(args.tabsNodeId);

		//var tabNode = dijit.byId(tabContainerNodeId);
		//tabsNode.addChild(tabNode);

		// SELECT THIS CHILD IN TABS
		//tabsNode.selectChild(tabNode);


        // LOAD SORIA AND FILEPICKER CSS
        this.loadCSS();

		// LOAD APPLICATIONS DND SOURCE AND TARGET INTO FRAMEWORK
		this.loadProjectTab();


	},


	//// GENERATE UNIQUE IDS FOR THE Projects PANE BY ADDING INCREMENTED PANE ID
	//setIds : function()
	//{
	//	// LAYOUT TEMPLATE
	//	var layout = this.layout;
	//	
	//	// SET UNIQUE IDS USING this.paneId
	//	if ( this.paneId.match(/^(.)(\D+)(\d+)$/) )
	//	{
	//		var firstLetter = this.paneId.match(/^(.)(\D+)(\d*)$/)[1];
	//		firstLetter = firstLetter.toUpperCase();
	//		var rest = this.paneId.match(/^(.)(\D+)(\d+)$/)[2];
	//		var number = this.paneId.match(/^(.)(\D+)(\d+)$/)[3];		
	//		layout.params.title = firstLetter + rest + " " + number;
	//	}
	//	
	//	layout.params.id = this.paneId + layout.params.id;
	//	for ( var i in layout.children )
	//	{
	//		layout.children[i].params.id = this.paneId + layout.children[i].params.id;
	//
	//		for ( var j in layout.children[i].children )
	//		{
	//			layout.children[i].children[j].params.id = this.paneId + layout.children[i].children[j].params.id;
	//		}
	//	}
	//
	//	this.layout = layout;		
	//},


	loadProjectTab: function ()
	{	

		var dialogId = dojo.dnd.getUniqueId();		

		// CREATE DIALOGUE
		var fileManagerDialog = new dojox.widget.Dialog ({
			id: dialogId,
			dimensions: [800,800],
			draggable: true,
			align: 'left',
			sizeDuration: 200,
			sizeMethod:       "combine",
			viewportPadding : "125",
			showTitle: "true",
			title: "File Manager"

		}, "fileManagerDialog" );		


		// SET CLASS FOR STYLE OF INSERTED PANE
		//dojo.attr(dojo.byId(dialogId), 'class', 'fileManager');

		// FADE OUT LOGIN WINDOW
		//dojo.fadeOut({ node: "loginDialogue", duration: 700 }).play();
		//setTimeout( "dijit.byId('loginDialogue').destroy()", 700);


		dojo.attr(fileManagerDialog.titleBar, 'class', 'fileManagerTitle');
		//dojo.marginBox(fileManagerDialog.titleBar).h = 30px;


		// getTitleHeight: function(){
		// summary: returns the height of the title dom node
		//  return dojo.marginBox(this.titleNode).h;    // Integer
		//},


		// center title bar
		//.dojoxDialogTitleBar {
		//	
		//}

		//dojo.attr(fileManagerDialog.titleBar, 'class', 'centered');
		//fileManagerDialog.titleBar.attr('class', 'fileManagerTitle');
		//fileManagerDialog.titleBar.attr("class", "centered");
		//fileManagerDialog.attr("class", "centered");




		var projects = this.projects();

		var sources = projects.sources;
		projects = projects.projects;


		// POPULATE PROJECTS PANE
		//for ( var i = 0; i < projects.length; i++ )
		for ( var i = 0; i < 1; i++ )
		{
			var project = projects[i].project;
			var owner = projects[i].owner;

			//var project = projectPath.match(/^\.*\/*(.+)$/)[1];



			//// CREATE 'PROJECTS' FILE DRAGGER
			//var paneNodeId = dojo.dnd.getUniqueId();
			//
			//
			////var paneNodeId = project + "-" + owner;
			//var paneNode = document.createElement('div');
			//paneNode.id = paneNodeId;
			//document.body.appendChild(paneNode);
			//
			//
			var containerNode = document.createElement('div');
			containerNode.id = paneNodeId + "-projectContainer";
			containerNode.innerHTML = project + "<BR>";
			//
			//
			//paneNode.appendChild(containerNode);




			var projectStore = new dojox.data.FileStore(
				{
					id: paneNodeId + "-fileStore",

					url: Agua.cgiUrl + "project.cgi?mode=fileSystem&sessionId=1228319084.3060.776&username=admin",
					pathAsQueryParam: true
				}
			);

			// SET FILE STORE path TO project
			projectStore.preamble = function()
			{
				this.store.path = this.arguments[0].path;                        
			};


			// GENERATE NEW FileDrag OBJECT
			var ProjectDrag = new plugins.files.FileSelector(
				{
					id: containerNodeId + "-ProjectDrag",
					style: "height: 270px; width: 100%;",
					store: projectStore,
					callback: this.callback
				}
			);


			// SET PATH FOR THIS PROJECT
			ProjectDrag.path = project;                    

			// START UP FileDrag
			ProjectDrag.startup();

			paneNode.appendChild(ProjectDrag.domNode);



		// SHOW LOGIN WINDOW
		// LATER: FIX THIS ERROR 'Error undefined running custom onLoad code'
		fileManagerDialog.setContent(ProjectDrag.domNode);
		fileManagerDialog.startup();
		fileManagerDialog.show();




			// APPEND FileDrag TO PANE
			//dojo.byId(sourceId).appendChild(paneNode);



		} // for loop on projects.length		


		// REMOVE EXISTING fileManagerDialog



		// SET fileManagerDialog
		this.fileManagerDialog = fileManagerDialog;

		//this.sourceId = sourceId;
		//this.targetId = targetId;
	}


	,


    loadCSS : function ()
    {        
        // THIS CSS FILE PROVIDES THE ICONS AND FORMATTING FOR INDIVIDUAL FILES/DIRECTORIES
	}


//	,
//	
//	
//    loadCSS : function ()
//    {        
//        // THIS CSS FILE PROVIDES THE ICONS AND FORMATTING FOR INDIVIDUAL FILES/DIRECTORIES
//		var cssFile1 = "dojo-1.5.0/dijit/themes/soria/soria.css";
//		var cssNode = document.createElement('link');
//		cssNode.type = 'text/css';
//		cssNode.rel = 'stylesheet';
//		cssNode.href = cssFile1;
//		cssNode.media = 'screen';
//		cssNode.title = 'loginCSS';
//		cssNode.id = "themeStyles";
//		document.getElementsByTagName("head")[0].appendChild(cssNode);
// 
//        // THIS CSS FILE PROVIDES THE ICONS AND FORMATTING FOR INDIVIDUAL FILES/DIRECTORIES
//		var cssFile3 = "plugins/files/FileDrag/FileDrag.css";
//		var cssNode = document.createElement('link');
//		cssNode.type = 'text/css';
//		cssNode.rel = 'stylesheet';
//		cssNode.href = cssFile3;
//		cssNode.media = 'screen';
//		cssNode.title = 'loginCSS';
//		cssNode.id = "widgetStyle";
//		document.getElementsByTagName("head")[0].appendChild(cssNode);
//
//
//        // THIS CSS FILE PROVIDES THE ICONS AND FORMATTING FOR INDIVIDUAL FILES/DIRECTORIES
//		var cssFile3 = "plugins/project/FileSelector/FileManager.css";
//		var cssNode = document.createElement('link');
//		cssNode.type = 'text/css';
//		cssNode.rel = 'stylesheet';
//		cssNode.href = cssFile3;
//		cssNode.media = 'screen';
//		cssNode.title = 'loginCSS';
//		cssNode.id = "widgetStyle";
//		document.getElementsByTagName("head")[0].appendChild(cssNode);
//
//    },
//
// 
// 

}); // end of Login

