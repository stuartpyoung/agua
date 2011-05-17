dojo.provide("plugins.report.Report");

/* CLASS SUMMARY: DISPLAY DATA TO USERS TO MODIFY OUTPUT IN REAL TIME
*/

// DnD
dojo.require("dojo.dnd.Source"); // Source & Target
dojo.require("dojo.dnd.Moveable");
dojo.require("dojo.dnd.Mover");
dojo.require("dojo.dnd.move");

// comboBox data store
dojo.require("dojo.data.ItemFileReadStore");
dojo.require("dijit.form.ComboBox");

// rightPane buttons
dojo.require("dijit.form.Button");

// PANE MODULES
dojo.require("plugins.core.BorderContainer");
dojo.require("plugins.core.ExpandoPane");

// PANE MODULES
dojo.require("dijit.layout.BorderContainer");
dojo.require("dojox.layout.ExpandoPane");
dojo.require("dijit.layout.ContentPane");

// INHERITS
dojo.require("plugins.core.Common");

// TEMPLATE MODULES
dojo.require("plugins.report.SNP");

dojo.declare( "plugins.report.Report", 
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	/////}}

// PANE ID 
paneId : null,
// PANE IDENTIFIER
//paneId : "reportPane",

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "report/templates/report.html"),

cssFiles : [
	"plugins/report/css/report.css",
	"dojo-1.5.0/dijit/themes/tundra/tundra.css",
	"dojox/layout/resources/ExpandoPane.css",
	"dojox/layout/resources/ExpandoPane.css",
	"dojox/layout/tests/_expando.css"
],

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// DEFAULT HEADINGS
tabs : [ "SNP" ],

// PANE WIDGETS
paneWidgets : null,

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

	// ADD ADMIN TAB TO TAB CONTAINER		
	Agua.tabs.addChild(this.mainTab);
	Agua.tabs.selectChild(this.mainTab);


//return;



	// CREATE HASH TO HOLD INSTANTIATED PANE WIDGETS
	this.paneWidgets = new Object;

	for ( var i = 0; i < this.tabs.length; i++ )
	{
		var tabPaneName = this.tabs[i];

		var thisObject = this;
		var attachNode = tabPaneName.toLowerCase() + "Pane";

		var tabPane = new plugins.report[tabPaneName](
			{
				filename: "454HCDiffs-headers-SNP.txt",
				project : this.project,
				workflow : this.workflow,
				attachNode : this[attachNode].domNode,
				gridAttachNode : this.gridAttachNode
			}
		);


		// REGISTER THE NEW TAB PANE IN this.paneWidgets 
		this.paneWidgets["plugins.admin." + tabPaneName] = tabPane;
	}


return;

},


//loadLeftPane: function ()
//{
//	
//	var reportOptions = new dojo.dnd.Source(
//		this.dragSource,
//		{
//			accept: ["none"],
//			copyOnly: true,
//		}
//	);
//	
//	if ( this.headings == null )
//	{
//		return;
//	}
//	reportOptions.insertNodes(false, this.headings);
//
//	// SET ONCLICK TO OPEN REPORT
//	var report = this;
//	var allNodes = reportOptions.getAllNodes();
//	for ( var i = 0; i < allNodes.length; i++ )
//	{
//		var node = allNodes[i];
//		var nodeClass = this.headings[i].type;
//		var applicationName = node.innerHTML;
//
//		// ADD CLASS FROM type TO NODE
//		dojo.addClass(node, nodeClass);
//		
//		// SET ONCLICK
//		node.onclick = function(e) {
//			
//			dojo.hitch( this, function(e)
//				{
//					this.loadRightPane(e);
//				}
//			);
//		}
//	}
//},
//

//loadMiddlePane : function()
//{
//	//
//	//// ADD REPORT CLASS
//	//dojo.addClass(this.middlePane.domNode, 'report');
//
//	// SET LEFT POSITION
//	//this.middlePane.domNode.style.left = this.middlePaneLeft;
//
//	// CREATE SNP REPORT
//	var className = name.replace(/\./g, '_');
//	var uniqueId = dijit.getUniqueId(className);
//	var snp = new plugins.report.SNP(
//		{
//			id: uniqueId,
//			filename: "Project1/Workflow1/454HCDiffs-headers-SNP.txt",
//			project : this.project,
//			workflow : this.workflow,
//			attachNode : this.middlePane.domNode,
//			gridAttachNode : this.rightPane.domNode
//		}
//	);
//
//	//middlePane.domNode.appendChild(snp.domNode);   
//},
//
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

	var reportObject = this;
	this.paneWidgets[paneName] = new plugins.report[paneName](
		{
			parentWidget: reportObject,
			tabContainer : reportObject.leftTabContainer
		}
	);
},
// GET THE LIST OF ADMIN PANES FROM THE SERVER BASED ON 
// AUTHENTICATION AND ADMIN STATUS OF USER 
getHeadings : function () {
	var url = Agua.cgiUrl + "report.cgi?";
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "getHeadings";

	var reportObject = this;

	// SEND TO SERVER
	dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			sync : true,
			handleAs: "json",
			putData: dojo.toJson(query),
			timeout: 15000,
			load: function(data)
			{
				headings = data;

				// SET this.headings
				reportObject.headings = headings;

			},
			error: function(response, ioArgs) {
				return response;
			}
		}
	);

	return headings;
}


}); // end of plugins.report.Report

