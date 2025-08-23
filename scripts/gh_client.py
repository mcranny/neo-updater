import json, base64, urllib.request, urllib.parse
API_GH = "https://api.github.com"

class GitHubClient:
    def __init__(self, token: str, repo: str, user_agent: str = "neo-snapshot/1.0"):
        self.token = token
        self.repo = repo
        self.ua = user_agent

    def _get(self, url, params=None, timeout=30):
        if params:
            url = f"{url}?{urllib.parse.urlencode(params)}"
        req = urllib.request.Request(url, headers={
            "Authorization": f"Bearer {self.token}",
            "Accept": "application/vnd.github+json",
            "User-Agent": self.ua
        }, method="GET")
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.getcode(), r.read()

    def _put(self, url, payload: dict, timeout=30):
        body = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(url, data=body, headers={
            "Authorization": f"Bearer {self.token}",
            "Accept": "application/vnd.github+json",
            "Content-Type": "application/json",
            "User-Agent": self.ua
        }, method="PUT")
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.getcode(), r.read()

    def get_contents(self, path: str, ref: str):
        url = f"{API_GH}/repos/{self.repo}/contents/{path}"
        code, body = self._get(url, params={"ref": ref})
        if code == 200:
            return json.loads(body.decode("utf-8"))
        if code == 404:
            return None
        raise RuntimeError(f"GET {url} -> {code}: {body[:200]}")

    def put_contents(self, path: str, message: str, content_bytes: bytes, sha: str|None, branch: str):
        content_b64 = base64.b64encode(content_bytes).decode("utf-8")
        payload = {"message": message, "content": content_b64, "branch": branch}
        if sha:
            payload["sha"] = sha
        url = f"{API_GH}/repos/{self.repo}/contents/{path}"
        code, body = self._put(url, payload)
        if code not in (200, 201):
            raise RuntimeError(f"PUT {url} -> {code}: {body[:400]}")
        return json.loads(body.decode("utf-8"))
