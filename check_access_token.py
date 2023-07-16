import tempfile
import os
import google_auth_oauthlib.flow

def get_access_token():
    """
    Retrieves the access token for Google API authentication.

    Returns:
        str: The access token for API authentication.

    Raises:
        IOError: If the current directory does not have write access to store the credentials.
    """

    print("Using access token. For info on how this is used see logging_info()")
    tf = tempfile.NamedTemporaryFile(delete=False).name
    try:
        with open(tf, "w"):
            pass
    except IOError:
        raise IOError("You are currently in a directory which doesn't have write access.\n"
                      "  In order to authenticate, we need to store the credentials in a file called '.httr-oauth'.\n"
                      "  Please set the working directory to a different directory where you have write access.")

    try:
        os.unlink(tf)
    except OSError:
        pass

    flow = google_auth_oauthlib.flow.InstalledAppFlow.from_client_secrets_file(
        "client_secrets.json", scopes=["https://www.googleapis.com/auth/userinfo.email"]
    )  # Replace "client_secrets.json" with your client secrets file

    credentials = flow.run_local_server(port=0)
    access_token = credentials.token

    return access_token

def check_access_token():
    """
    Checks if the access token file exists and retrieves the access token.

    Returns:
        str or None: The access token if it exists, otherwise None.

    Example:
        token = check_access_token()
        if token:
            print("Access token:", token)
        else:
            print("No access token found.")
    """

    if os.path.exists("ieugwasr_oauth"):
        return get_access_token()
    else:
        return None
