# advancedcg_report
Repository for the lecture course "Advanced Computer Graphics" in University of Tsukuba

## トラブルシューティング
 - Visual Studioはgithub連携機能がありますが，このリポジトリには複数のVisual Studioソリューションとソースコード以外のファイルが多く含まれているのでうまくいかない可能性があります．別途リポジトリをクローンするかDownload ZIPで全体をダウンロード・解凍後，srcフォルダの各フォルダ内に含まれる*.slnファイルをVisual Studioで開くようにしてください．

 - Visual Studioで「error MSB8036: Windows SDK バージョン 10.0.xxxxxx.0 が見つかりませんでした。」というエラーが出てビルドできない．  
　⇒ 「プロジェクト」メニューから「プロジェクトの再ターゲット」で「Windows SDK バージョン:」 のところにその環境で対応するバージョンが出るので，問題なければそのまま「OK」をクリック
 - Visual Studioで実行しても特定のプロジェクトしか実行されない場合は，「プロジェクト」メニューから「スタートアップ プロジェクトの構成」で「現在の選択」にチェックしてください．

 - 自身のPC環境(Windows)で動かしたい人向けに講義資料フォルダ(SharePointの方)に「Visual Studio Community インストール方法」というPDFを用意してあるので参考にしてください．
