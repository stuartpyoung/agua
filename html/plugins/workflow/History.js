dojo.provide("plugins.workflow.History");

// ALLOW THE USER TO ADD, REMOVE AND MODIFY PARAMETERS

dojo.require("dijit.layout.ContentPane");

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dojo.parser");
dojo.require("plugins.core.Common");

// HAS A
//dojo.require("plugins.workflow.HistoryRow");

dojo.declare(
    "plugins.workflow.History",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/history.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ dojo.moduleUrl("plugins") + "/workflow/css/history.css" ],

// PARENT WIDGET
parentWidget : null,

// ARRAY OF CHILD WIDGETS
childWidgets : null,

// isVALID BOOLEAN: ALL PARAMETERS ARE VALID
isValid : null,

// CORE WORKFLOW OBJECTS
core : null,

/////}

constructor : function(args)
{
	this.core = args.core;

	// LOAD CSS
	this.loadCSS();		
},

postCreate : function()
{

	this.startup();
},


// DO inherited, LOAD ARGUMENTS AND ATTACH THE MAIN TAB TO THE ATTACH NODE
startup : function ()
{

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 


	// ADD TO TAB CONTAINER		
	this.attachNode.addChild(this.mainTab);
	this.attachNode.selectChild(this.mainTab);

	// GET WORKFLOW HISTORY AND DISPLAY IN HISTORY TAB
	this.showHistory();
},

// GET WORKFLOW HISTORY AND DISPLAY IN HISTORY TAB
showHistory : function ()
{

	// SELECT HISTORY TAB IN TAB CONTAINER
	//this.controlPane.selectChild(this.history);

	// EMPTY CURRENT TABLE IF PRESENT
	while ( this.historyTable.firstChild )
	{
		this.historyTable.removeChild(this.historyTable.firstChild);
	}

	// GET URL 
	var url = Agua.cgiUrl + "workflow.cgi";

	// GENERATE QUERY JSON FOR THIS WORKFLOW IN THIS PROJECT
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "getHistory";

	var thisObject = this;
	dojo.xhrPut(
		{
			url: url,
			putData: dojo.toJson(query),
			handleAs: "json",
				//handleAs: "json-comment-optional",
			sync: false,
			handle: function(response)
			{
				// BUILD TABLE
				var table = document.createElement('table');
				thisObject.historyTable.appendChild(table);
				dojo.addClass(table, 'historyTable');

				// SET THE NODE CLASSES BASED ON STATUS
				for ( var i = 0; i < response.length; i++ )
				{
					var tr = document.createElement('tr');
					table.appendChild(tr);

					var td = document.createElement('td');
					tr.appendChild(td);

					var project = response[i][0].project;
					var workflow = response[i][0].workflow;
					var historyPane = new plugins.workflow.HistoryPane(
						{
							project: project,
							workflow : workflow,
							rows: response[i]
						}
					);
					td.appendChild(historyPane.domNode); 
				}
			}
		}
	);
},



downloadHistory: function (event)
{
	// STOP EVENT
	event.stopPropagation();

	// GENERATE HISTORY QUERY AND SEND TO SERVER
	var query = "?username=" + Agua.cookie('username');
	query += "&sessionId=" + Agua.cookie('sessionId');
	query += "&mode=downloadHistory";
	var url = Agua.cgiUrl + "download.cgi";

	// DO AN IFrame REQUEST TO DOWNLOAD THE FILE
	/* NB: dojo.io.iframe.send SUPPORTS ONLY 'GET' OR 'POST */
	var args = {
		method: "POST",
		url: url + query,
		//content: {},
		handleAs: "html",
		timeout: 10000
		//load: dojo.hitch(this, "onDownloadComplete"),
		//error: dojo.hitch(this, "onDownloadError")
	};
	dojo.io.iframe.send(args);
}



}); // plugins.workflow.History

