dojo.provide("plugins.core.Updater");

// UPDATE SUBSCRIBING OBJECTS, E.G., IN RESPONSE TO INFORMATION
// CHANGES THAT AFFECT THE SUBSCRIBING OBJECTS
dojo.declare( "plugins.core.Updater",
    null,
{    
// HASH OF LOADED CSS FILES
loadedCssFiles : null,

updates : null,

constructor : function () {
	this.startup();
},

startup : function () {
	this.updates = new Object;
},

subscribe : function (subscriber, subscription) {

	// CHECK SUBSCRIBER IMPLEMENTS SUBSCRIPTION METHOD
	if ( subscriber[subscription] == null )
	{
		console.error("core.Updater.subscribe    subscriber " + subscriber + " doesn't implement subscription method: " + subscription);
		return;
	}

	if ( this.updates[subscription] == null )
		this.updates[subscription] = new Array;

	for ( currentSubscriber in this.updates[subscription] )
	{
		if ( currentSubscriber == subscriber )
		{
			return;	
		}
	}

	this.updates[subscription].push(subscriber);

},

update : function (subscription, args) {

	var subscribers = this.getSubscribers(subscription);
	for ( var i in subscribers )
	{
		var subscriber = subscribers[i];
		subscriber[subscription](args);
	}
},

getSubscribers : function (subscription) {
	return this.updates[subscription]
},

unsubscribe : function (subscriber, subscription) {

	// CHECK SUBSCRIBER IMPLEMENTS SUBSCRIPTION METHOD
	if ( subscriber[subscription] == null )
	{
		console.error("core.Updater.unsubscribe    subscriber " + subscriber + " doesn't implement subscription method: " + subscription);
		return;
	}

	if ( this.updates[subscription] == null )
		this.updates = new Array;

	for ( var i = 0; i < this.updates[subscription].length; i++ )
	{
		if ( this.updates[subscription][i] == subscriber )
		{
			this.updates[subscription].splice(i, 1);
			break;	
		}
	}
}



});

