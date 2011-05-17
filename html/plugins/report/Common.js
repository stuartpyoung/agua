
dojo.provide("plugins.report.Common");

// USEFUL METHODS FOR REPORT SUBCLASSES

dojo.declare( "plugins.report.Common", null,
{

 	// CONSTRUCTOR	
	constructor : function(args) {

	},

	// saveReport
	//
	// SEND ALL VALUES OF ALL ELEMENTS IN THE TEMPLATE TO THE SERVER
	//
	// Send report details to server with POST
	//// Get a reference to the Editor instance (through dijit.byId(), a jsid,
	//// or whatnot), then call editor.getValue(), which should give you the
	//// HTML markup it used.   That can then be post/put through dojo.xhr back
	//// to a service to save as HTML however you please.   As for making sure
	//// it is utf-8 encoded, set the post/put ContentEncoding to utf 8 and
	//// make sure your backend service saves the file in UTF-8
	//// http://en.wikipedia.org/wiki/UTF-8
	saveReport : function()
	{

		// STORE JSON QUERY
		var json = this.getData();

		// ADD mode AND class
		var query = new Object;
		query.mode = "saveReport";
		query["class"] = "Report::SNP";   // ['class'] avoids use of predefined 'class' operator

		// ADD USER AUTHENTICATION INFO AND REPORT INFO
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.project = Agua.cookie('project');
		query.workflow = Agua.cookie('workflow');
		query.report = Agua.cookie('report');
		query.json = json;


		// POST DATA TO SERVER
		dojo.xhrPut(
			{
				url: Agua.cgiUrl + "report.cgi",
				contentType: "text",
				//handleAs: "json",
				//handleAs: "json-comment-filtered",
				putData: dojo.toJson(query),
				timeout: 20000,
				load: function(response, ioArgs) {
					return response;
				},
				error: function(response, ioArgs) {
					return response;
				}
			}
		);

	},


	// getData
	// 
	// GET VALUES OF ALL ELEMENTS IN THE TEMPLATE
	//
	getData : function ()
	{
		var data = new Object;
		for ( var name in this.elementObjects )
		{
			if ( ! this.elementObjects[name].valueFunction )
			{
				data[name] = '';
			}
			else
			{
				var value = this.elementObjects[name].valueFunction();
				if ( ! value )
				{
					value = false;
				}
				data[name] = value;
			}
		}

		return data;
	}







}); // plugins.report.Report



