dojo.provide("plugins.core.form.Template");

// REFERENCE SITES

// Object Oriented concepts and inheritance
// http://dojotoolkit.org/book/dojo-book-0-4/part-3-dojo-programming-model/object-oriented-concepts-and-inheritance

// BIOMART
//http://www.ensembl.org/biomart/martview/e61a39e4e8e306354f8f2a7a70bbc53c/e61a39e4e8e306354f8f2a7a70bbc53c

// DOJO DEMOS
//http://localhost:8080/Bioptic0.2.5/html/dojo-1.5.0/demos/
//http://localhost:8080/Bioptic0.2.5/html/dojo-1.5.0/dojox/form/tests/test_SelectStack.html

//Real World Dojo part 5: Custom Components
//http://www.dotnetmafia.com/blogs/jamesashley/archive/2008/10/28/761.aspx

dojo.declare(
    "plugins.core.form.Template",
    [dijit._Widget, dijit._Templated],
    {
//        //Path to the template of this widget. dijit._Templated uses this to
//        //snatch the template from the named file via a synchronous call.
//        // This provides the relative path to the template file for the widget. Note that
//        // fetching the template for a widget requires a synchronous network call, although
//        // Dojo will cache the template string after it is fetched (unless a custom build
//        // with interned template strings is used).
//        templatePath: null,
//
//        // Calls dijit._Templated.widgetsInTemplate
//        // If widgets are defined inside of the template (path or string), this
//        // value should be specified and set to true so that the Dojo parser will
//        // find and instantiate those widgets. This value is false by default.
//        // While this property is easy to overlook, know that defining widgets
//        // inside of other widget?s templates can often be very useful.
//        widgetsInTemplate : true,
//        
//        // SLOTS TO BE POPULATED BY THE DESCENDENT OF THIS CLASS
//        title : '',
//        id : '',
//		filename: '',
//        url : '',       // url to save report to
//        widgets : [],
//        widgetObjects : [],
//        cssFiles : [],
//
//        //Any initialization code would go here in the constructor. dijit._Widget and
//        //dijit._Templated do not have parameters in their constructors, so
//        //there wouldn't be any multiple-inheritance complications
//        //if you were to include some paramters here.
//        constructor : function (args)
//        {
//
//
//    
//            this.id         = args.id;
//			this.filename   = args.filename;
//			this.url        = args.url;
//            
//            //// LOAD CSS FOR WIDGETS
//            //this.loadCSS();
//            
//            // SET UNIQUE IDS FOR EACH WIDGET
//            this.setWidgetIds();
//        },
//        
//
//        //Inherited from dijit._Widget and called just before template
//        //instantiation in buildRendering. This method is especially useful
//        //for manipulating the template before it becomes visible.
//        postMixInProperties: function()
//        {
//        },
//
//
//        //You can override this method to manipulate widget once it is
//        //placed in the UI, but be warned that any child widgets contained
//        //in it may not be ready yet.        
//        postCreate: function() {},
//
//
//        
//        //loadCSS : function()
//        //{
//        //    
//        //    // LOAD CSS
//        //    var cssFiles = this.cssFiles;
//        //    
//        //    for ( var i in cssFiles )
//        //    {
//        //        var cssFile = cssFiles[i];
//        //        
//        //        var cssNode = document.createElement('link');
//        //        cssNode.type = 'text/css';
//        //        cssNode.rel = 'stylesheet';
//        //        cssNode.href = cssFile;
//        //        cssNode.media = 'screen';
//        //        document.getElementsByTagName("head")[0].appendChild(cssNode);
//        //    }
//        //},
//        
//
//        setWidgetIds : function ()
//        {
//
//            var uniqueId = 0;
//            for (var index in this.widgets)
//            {
//                // GET THE NAME OF THE WIDGET
//				var name = this.widgets[index];
//                
//				// GET UNIQUE ID FOR WIDGET
//				var uniqueId = dijit.getUniqueId(this.id);
//                
//				// SET UNIQUE ID FOR WIDGET
//				this[name] = uniqueId;
//
//                // ADD THIS WIDGET'S name AND id TO A widgetObject OBJECTS IN THE
//                // this.widgetObjects HASH, USING name AS THE KEY
//                var widgetObject = new Object;
//                widgetObject.name = name;
//                widgetObject.id = uniqueId;
//                this.widgetObjects[name] = widgetObject;
//            }
//			
//
//        },
//
//
//
//		// NB: ONLY RETURNS A FIRST LAYER CHILD (NOT RECURSIVE ON CHILDREN OF CHILDREN)
//		getElementById : function (id)
//		{
//			
//			var childNodes = this.domNode.childNodes;
//			for ( var i = 0; i < childNodes; i++ )
//			{
//				if ( childNodes[i].id == id )
//				{
//					return childNodes[i];
//				}
//			}
//
//			return null;
//		},
//	    
//        //Called after the widget's children and all other widgets on the
//        //page have been created. Provides an opportunity to manipulate child
//        //widgets before they're displayed.
//        // NOW THAT ALL THE CHILD WIDGETS HAVE BEEN CREATED, ADD EACH CHILD WIDGET
//        // TO ITS RELATED widgetObject IN THE GLOBAL this.widgetObjects HASH
//        
//		// TO DO: FIX NON-WIDGET CASE BY RETURNING THE DOM NODE
//		startup : function ()
//        {
//
//            for ( var i in this.widgetObjects )
//            {
//                var widgetObject = this.widgetObjects[i];
//                var name = widgetObject.name;
//                var id = widgetObject.id;
//                
//                // ADD WIDGET TO widgetObject
//                var widget = dijit.byId(id);
//
//				// *** FIX THIS LATER ***
//				if ( ! widget )
//				{
//					var node = this.getElementById(id);
//					widget = node;
//				}
//                widgetObject.widget = widget;
//                this.widgetObjects[i] = widgetObject;
//
//                // RUN SPECIFIC FUNCTION FOR WIDGET TO ADD SPECIFIC PROPERTIES.
//                // USES THE FACT THAT this[functionName] CAN BE A FUNCTION
//                var widgetFunction = name + "Function";
//                this[widgetFunction](widget, name, id);
//            }
//        },
//        
//
//		getWidget : function (name)
//		{
//
//            //for ( var i = 0; i < this.widgetObjects.length; i++ )
//            //{
//            //    
//            //}
//
//			if ( ! this.widgetObjects )
//			{
//				return;
//			}
//			
//			//console.dir(this.widgetObjects[name]);
//
//			var widgetObject = this.widgetObjects[name];
//			if ( ! widgetObject )
//			{
//				return;
//			}
//
//            //return widgetObject.widget ? widgetObject.widget : null;
//			return widgetObject.widget;
//		},
//		
//
//		getNode : function (name)
//		{
//
//			if ( ! this.widgetObjects )
//			{
//				return;
//			}
//			
//			var widgetObject = this.widgetObjects[name];
//            var id = widgetObject.id;
//            if ( ! id )
//            {
//            }
//
//            return dojo.byId(id);
//		},
//		
//        
//        
//        loadWidgetValues : function (response, ioArgs)
//        {
//
//			for ( var name in response )
//			{
//                var value = response[name];
//				this.widgetObjects[name].value(this.widgetObjects[name], value);
//				//this.widgetObjects[name].value(this.widgetObjects[name], "'" + value + "'");
//			}
//        },
//
//
//        getValue : function (name)
//        {
//            
//            if ( ! this.widgetObjects[name] )
//            {
//                return '';
//            }
//            else if ( ! this.widgetObjects[name].value )
//            {
//                return '';
//            }
//
//            return this.widgetObjects[name].value(this.widgetObjects[name]);
//        },
//
//
//        setValue : function (name, value)
//        {
//            
//            if ( ! this.widgetObjects[name] )
//            {
//                return;
//            }
//            this.widgetObjects[name].value(this.widgetObjects[name], value);
//        },
//        
//
//        // EXAMPLE OF xhrPost
//        // http://www.dojoforum.com/2007/10/11/dojo-example-xhrget-and-xhrpost
//
//		loadReport : function()
//		{
//            	
//			// GENERATE QUERY
//			var query = new Object;
//			query.username = Agua.cookie('username');
//			query.sessionId = Agua.cookie('sessionId');
//			query.project = Agua.cookie('project');
//			query.workflow = Agua.cookie('workflow');
//			query.report = Agua.cookie('Report1');
//			query.mode = "loadReport";
//			query.module = "Report::SNP";
//	
//			// POST DATA TO SERVER
//            var template = this;
//			dojo.xhrPost(
//				{
//					//url: this.url,
//                    url : "../plugins/report/templates/SNP.json",
//					contentType: "text",
//					handleAs: "json",
//					postData: dojo.toJson(query),
//					//handleAs: "json-comment-filtered",
//					timeout: 2000,
//					load: function(response, ioArgs) {
//                        template.loadWidgetValues(response, ioArgs);
//                    },
//					error: function(response, ioArgs) {
//						return response;
//					}
//				}
//			);
//			
//		},
//
//
//        getData : function ()
//        {
//            var data = new Object;
//			for ( var name in this.widgetObjects )
//			{
//                if ( ! this.widgetObjects[name].value )
//                {
//                    data[name] = '';
//                }
//                else
//                {
//                    var value = this.widgetObjects[name].value(this.widgetObjects[name]);
//                    if ( ! value )
//                    {
//                        value = false;
//                    }
//                    data[name] = value;
//                }
//			}
//    
//            return data;
//        },
//        
//		
//		// Send report details to server with POST
//        //// Get a reference to the Editor instance (through dijit.byId(), a jsid,
//        //// or whatnot), then call editor.getValue(), which should give you the
//        //// HTML markup it used.   That can then be post/put through dojo.xhr back
//        //// to a service to save as HTML however you please.   As for making sure
//        //// it is utf-8 encoded, set the post/put ContentEncoding to utf 8 and
//        //// make sure your backend service saves the file in UTF-8
//        //// http://en.wikipedia.org/wiki/UTF-8
//		saveReport : function()
//		{
//            
//            // STORE JSON QUERY
//            var query = this.getData();
//
//			// ADD mode AND class
//			query["mode"] = "saveReport";
//			query["class"] = "Report::SNP";
//	
//            // ADD USER AUTHENTICATION INFO AND REPORT INFO
//			query.username = Agua.cookie('username');
//			query.sessionId = Agua.cookie('sessionId');
//			query.project = Agua.cookie('project');
//			query.workflow = Agua.cookie('workflow');
//			query.report = Agua.cookie('report');
//            
//	
//			// POST DATA TO SERVER
//			dojo.rawXhrPost(
//				{
//					url: this.url,
//					contentType: "text",
//					//handleAs: "json",
//					//handleAs: "json-comment-filtered",
//					postData: dojo.toJson(query),
//					timeout: 2000,
//					load: function(response, ioArgs) {
//						return response;
//					},
//					error: function(response, ioArgs) {
//						return response;
//					}
//				}
//			);
//			
//		}
    }
); // plugins.core.form.Template
